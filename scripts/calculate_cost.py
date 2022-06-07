#!/usr/bin/env python3
import argparse
import gzip
import json
import sys
import urllib.error
import urllib.parse
import urllib.request
import re
from io import BytesIO

# pip install python-dateutil
import dateutil.parser

GCE_MACHINE_TYPES_URL = "http://cloudpricingcalculator.appspot.com/static/data/pricelist.json"
TOTAL_WORKFLOW_COST = 0
TOTAL_RUN_HOURS = 0
CUSTOM_MACHINE_CPU = "CP-COMPUTEENGINE-CUSTOM-VM-CORE"
CUSTOM_MACHINE_RAM = "CP-COMPUTEENGINE-CUSTOM-VM-RAM"
CUSTOM_MACHINE_EXTENDED_RAM = "CP-COMPUTEENGINE-CUSTOM-VM-EXTENDED-RAM"
CUSTOM_MACHINE_CPU_PREEMPTIBLE = "CP-COMPUTEENGINE-CUSTOM-VM-CORE-PREEMPTIBLE"
CUSTOM_MACHINE_RAM_PREEMPTIBLE = "CP-COMPUTEENGINE-CUSTOM-VM-RAM-PREEMPTIBLE"
CUSTOM_MACHINE_EXTENDED_RAM_PREEMPTIBLE = "CP-COMPUTEENGINE-CUSTOM-VM-EXTENDED-RAM-PREEMPTIBLE"
CUSTOM_MACHINE_TYPES = [CUSTOM_MACHINE_CPU,
                        CUSTOM_MACHINE_RAM,
                        CUSTOM_MACHINE_EXTENDED_RAM,
                        CUSTOM_MACHINE_CPU_PREEMPTIBLE,
                        CUSTOM_MACHINE_RAM_PREEMPTIBLE,
                        CUSTOM_MACHINE_EXTENDED_RAM_PREEMPTIBLE]


# load the US pricing for both persistent disk and compute engine
def get_gce_pricing():
    response = urllib.request.urlopen(GCE_MACHINE_TYPES_URL)
    data = response.read()

    if response.info().get('Content-Encoding') == 'gzip':
        buf = BytesIO(data)
        f = gzip.GzipFile(fileobj=buf)
        data = f.read()

    pricing = json.loads(data)

    data = {}

    for k, v in pricing.items():
        if k == "gcp_price_list":
            for k2, v2 in v.items():
                if k2.startswith("CP-COMPUTEENGINE-VMIMAGE"):
                    data[k2.replace("CP-COMPUTEENGINE-VMIMAGE-", "").lower()] = get_us_keys(v2)
                if k2.startswith("CP-COMPUTEENGINE-STORAGE-PD"):
                    data[k2] = get_us_keys(v2)
                if k2.startswith("CP-COMPUTEENGINE-LOCAL-SSD"):
                    data[k2] = get_us_keys(v2)
                if k2.startswith("CP-COMPUTEENGINE-LOCAL-SSD-PREEMPTIBLE"):
                    data[k2] = get_us_keys(v2)
                if k2.startswith("CP-COMPUTEENGINE-STORAGE-PD-CAPACITY"):
                    data[k2] = get_us_keys(v2)
                if k2 in CUSTOM_MACHINE_TYPES:
                    data[k2] = get_us_keys(v2)

    return data

def get_us_keys(dict):
    us_cost = count = 0
    for k, v in dict.items():
        if k.startswith("us"):
            count += 1
            us_cost += v
    if count == 0:
        print("The dictionary doesn't contain a US cost object")
        exit(0)
    return us_cost / count

def extract_machine_type(call_info):
    # First, look in the executionEvents for the type of machine allocated:
    if 'executionEvents' in call_info:
        for event in call_info['executionEvents']:
            description = event['description']
            match = re.search("^Worker.*assigned in \"(\S+)\" on a \"(\S+)\" machine", description)
            if match:
                return match.groups()[1]

    # If here, the machine type was not found in the executionEvents, look in jes
    if 'jes' in call_info and 'machineType' in call_info['jes']:
        full_machine = call_info['jes']['machineType']
        if full_machine.startswith("custom"):
            return "custom"
        elif full_machine.startswith("predefined"):
            return "custom"
        else:
            return full_machine.split("/")[1]
    else:
        return "unknown"


def get_disk_info(metadata):
    if "runtimeAttributes" in metadata and "disks" in metadata['runtimeAttributes']:
        bootDiskSizeGb = 0.0
        if "bootDiskSizeGb" in metadata['runtimeAttributes']:
            bootDiskSizeGb = float(metadata['runtimeAttributes']['bootDiskSizeGb'])
        # Note - am lumping boot disk in with requested disk.  Assuming boot disk is same type as requested.
        # i.e. is it possible that boot disk is HDD when requested is SDD.
        (name, disk_size, disk_type) = metadata['runtimeAttributes']["disks"].split()
        return {"size": float(disk_size) + bootDiskSizeGb, "type": "PERSISTENT_" + disk_type}
    else:
        # we can't tell disk size in this case so just return nothing
        return {"size": float(0), "type": "PERSISTENT_SSD"}


def was_preemptible_vm(metadata):
    if "runtimeAttributes" in metadata and "preemptible" in metadata['runtimeAttributes']:
        pe_count = int(metadata['runtimeAttributes']["preemptible"])
        attempt = int(metadata['attempt'])

        return attempt <= pe_count
    else:
        # we can't tell (older metadata) so conservatively return false
        return False

def used_cached_results(metadata):
    return "callCaching" in metadata and metadata["callCaching"]["hit"]

def calculate_runtime(call_info, ignore_preempted):
    # get start (start time of VM start) & end time (end time of 'ok') according to metadata
    start = None
    end = None

    # give a runtime of 0 for preempted jobs so they have no cost associated with them
    if was_preempted(call_info) and ignore_preempted:
        return 0

    if 'executionEvents' in call_info:
        for x in call_info['executionEvents']:
            y = x['description']

            if call_info['backend'] == 'PAPIv2':
                if y.startswith("PreparingJob"):
                    start = dateutil.parser.parse(x['startTime'])
                if y.startswith("Worker released"):
                    end = dateutil.parser.parse(x['endTime'])
            else:
                if y.startswith("start"):
                    start = dateutil.parser.parse(x['startTime'])
                if y.startswith("ok"):
                    end = dateutil.parser.parse(x['endTime'])

    # if we are preempted or if cromwell used previously cached results, we don't even get a start time from JES.
    # if cromwell was restarted, the start time from JES might not have been written to the metadata.
    # in either case, use the Cromwell start time which is earlier but not wrong.
    if start is None:
        start = dateutil.parser.parse(call_info['start'])

    # if we are preempted or if cromwell used previously cached results, we don't get an endtime from JES right now.
    # if cromwell was restarted, the start time from JES might not have been written to the metadata.
    # in either case, use the Cromwell end time which is later but not wrong
    if end is None:
        end = dateutil.parser.parse(call_info['end'])

    # The minimum runtime is 1 minute, after that it's by the second.
    # so if the task ran for 30 seconds, you pay for 1 minute.  If it ran for 1:01 minute, you pay for 1:01 minute
    elapsed = end - start
    seconds = elapsed.days * 24 * 60 * 60 + elapsed.seconds
    run_seconds = max(60.0, seconds)
    return run_seconds

def was_preempted(call_info):
    # We treat Preempted and RetryableFailure the same.  The latter is a general case of the former
    return call_info['executionStatus'] in ['Preempted', 'RetryableFailure']


def calculate_cost(metadata, ignore_preempted, only_total_cost, print_header):
    # set up pricing information
    pricing = get_gce_pricing()
    ssd_cost_per_gb_per_month = float(pricing["CP-COMPUTEENGINE-STORAGE-PD-SSD"])
    ssd_cost_per_gb_hour = (ssd_cost_per_gb_per_month / (24 * 365 / 12))

    local_ssd_cost_per_gb_per_month = float(pricing["CP-COMPUTEENGINE-LOCAL-SSD"])
    local_ssd_cost_per_gb_hour = (ssd_cost_per_gb_per_month / (24 * 365 / 12))

    pe_local_ssd_cost_per_gb_per_month = float(pricing["CP-COMPUTEENGINE-LOCAL-SSD-PREEMPTIBLE"])
    pe_local_ssd_cost_per_gb_hour = (ssd_cost_per_gb_per_month / (24 * 365 / 12))

    hdd_cost_per_gb_per_month = float(pricing["CP-COMPUTEENGINE-STORAGE-PD-CAPACITY"])
    hdd_cost_per_gb_hour = (hdd_cost_per_gb_per_month / (24 * 365 / 12))

    disk_costs = {"PERSISTENT_SSD"      : ssd_cost_per_gb_hour,
                  "PERSISTENT_HDD"      : hdd_cost_per_gb_hour,
                  "PERSISTENT_LOCAL"    : local_ssd_cost_per_gb_hour,
                  "PE_PERSISTENT_LOCAL" : pe_local_ssd_cost_per_gb_hour,
                  }

    if print_header and not only_total_cost:
        # print out a header
        print("\t".join(
            ["task_name", "status", "machine_type", "cpus", "mem_gbs",
             "total_hours", "cpu_cost_per_hour", "cpu_cost", "mem_cost_per_hour", "mem_cost",
             "pe_total_hours", "pe_cpu_cost_per_hour", "pe_cpu_cost", "pe_mem_cost_per_hour", "pe_mem_cost",
             "failed_pe_total_hours", "failed_pe_cpu_cost", "failed_pe_mem_cost",
             "disk_type", "disk_size", "disk_gb_hours", "disk_cost",
             "failed_pe_ssd_gb_hours", "failed_pe_ssd_cost",
             "total_cost"]))

    # iterate through the metadata file for each call
    for k, v in metadata['calls'].items():
        task_name = k

        total_hours = 0
        pe_total_hours = 0
        failed_pe_total_hours = 0
        cpus = 0
        mem_gbs = 0
        machine_type = "unknown"
        complete = True
        disk_info = get_disk_info({})

        for call_info in v:
            # this is a subworkflow, recursively calculate cost on workflow metadata
            if 'subWorkflowMetadata' in call_info:
                calculate_cost(call_info['subWorkflowMetadata'], ignore_preempted, only_total_cost, False)
            else:
                # only process things that are not in flight
                if call_info['executionStatus'] in ['Running', 'NotStarted', 'Starting']:
                    complete = False
                else:
                    if call_info['executionStatus'] in ['Failed']:
                        complete = False

                    if machine_type == "unknown":
                        machine_type = extract_machine_type(call_info)

                    pe_vm = was_preemptible_vm(call_info)
                    disk_info = get_disk_info(call_info)

                    run_hours = calculate_runtime(call_info, ignore_preempted) / (60.0 * 60.0)

                    # for preemptible VMs, separately tally successful tasks vs ones that were preempted
                    if pe_vm:
                        if was_preempted(call_info):
                            # If Compute Engine terminates a preemptible instance less than 10 minutes after it is created,
                            # you are not billed for the use of that virtual machine instance
                            if run_hours < (10.0 / 60.0):
                                run_hours = 0
                            failed_pe_total_hours += run_hours
                        else:
                            pe_total_hours += run_hours
                    else:
                        total_hours += run_hours

        # Runtime parameters are the same across all calls; just pull the info from the first one
        if 'runtimeAttributes' in v[0]:
            if 'cpu' in v[0]['runtimeAttributes']:
                cpus += int(v[0]['runtimeAttributes']['cpu'])
            if 'memory' in v[0]['runtimeAttributes']:
                mem_str = v[0]['runtimeAttributes']['memory']
                mem_gbs += float(mem_str[:mem_str.index(" ")])

        if complete:
            status = "complete"
        else:
            status = "incomplete"

        if machine_type != "custom" and machine_type not in pricing:
            if "n2" in machine_type:
                machine_type = machine_type.replace("n2", "n1")
                if machine_type not in pricing:
                    machine_type = "unknown"
            else:
                machine_type = "unknown"

        if machine_type == "unknown":
            cpu_cost_per_hour = 0
            pe_cpu_cost_per_hour = 0
            mem_cost_per_hour = 0
            pe_mem_cost_per_hour = 0
        elif machine_type == "custom":
            cpu_cost_per_hour = pricing[CUSTOM_MACHINE_CPU] * cpus
            pe_cpu_cost_per_hour = pricing[CUSTOM_MACHINE_CPU_PREEMPTIBLE] * cpus
            mem_cost_per_hour = pricing[CUSTOM_MACHINE_RAM] * mem_gbs
            pe_mem_cost_per_hour = pricing[CUSTOM_MACHINE_RAM_PREEMPTIBLE] * mem_gbs
        else:
            cpu_cost_per_hour = pricing[machine_type]
            pe_cpu_cost_per_hour = pricing[machine_type + "-preemptible"]
            mem_cost_per_hour = 0
            pe_mem_cost_per_hour = 0

        cpu_cost = total_hours * cpu_cost_per_hour
        failed_pe_cpu_cost = failed_pe_total_hours * pe_cpu_cost_per_hour
        pe_cpu_cost = pe_total_hours * pe_cpu_cost_per_hour

        #
        # NOTE -- local ssds have a different price when used in preemptible VMs.  However, to implement this all the disk calculations
        # need to be moved from the task level (where it is now) to the call level since each call could be preemptible or not
        # Then we can decide to use PERSISTENT_LOCAL or PE_PERSISTENT_LOCAL
        #
        disk_cost_per_gb_hour = disk_costs[disk_info["type"]]

        disk_gb_hours = disk_info["size"] * (total_hours + pe_total_hours)
        disk_cost = disk_gb_hours * disk_cost_per_gb_hour

        failed_pe_disk_gb_hours = disk_info["size"] * failed_pe_total_hours
        failed_pe_disk_cost = failed_pe_disk_gb_hours * disk_cost_per_gb_hour

        mem_cost = total_hours * mem_cost_per_hour
        pe_mem_cost = pe_total_hours * pe_mem_cost_per_hour
        failed_pe_mem_cost = failed_pe_total_hours * pe_mem_cost_per_hour

        total_cost = cpu_cost + \
                     pe_cpu_cost + \
                     failed_pe_cpu_cost + \
                     disk_cost + \
                     failed_pe_disk_cost + \
                     mem_cost + \
                     pe_mem_cost + \
                     failed_pe_mem_cost

        # accumalate total workflow cost
        global TOTAL_WORKFLOW_COST
        TOTAL_WORKFLOW_COST += total_cost

        global TOTAL_RUN_HOURS
        TOTAL_RUN_HOURS += total_hours + pe_total_hours

        if not only_total_cost:
            out = (
                task_name, status, machine_type, cpus, mem_gbs,
                total_hours, cpu_cost_per_hour, cpu_cost, mem_cost_per_hour, mem_cost,
                pe_total_hours, pe_cpu_cost_per_hour, pe_cpu_cost, pe_mem_cost_per_hour, pe_mem_cost,
                failed_pe_total_hours, failed_pe_cpu_cost, failed_pe_mem_cost,
                disk_info["type"], disk_info["size"], disk_gb_hours, disk_cost,
                failed_pe_disk_gb_hours, failed_pe_disk_cost,
                total_cost)
            print('\t'.join(map(str, out)))


def compare(old, new):
    """Fail when NEW total exceeds OLD total by > 5%."""
    def total(cost_file):
        with open(cost_file, 'r') as input:
            lines = input.readlines()
        for line in lines:
            fields = line.split()
            if len(fields) == 3 and fields[0] == 'Total' and fields[1] == 'Cost:':
                return int(float(fields[2]) * 10000) / 10000.0
        return None
    old_cost = total(old)
    new_cost = total(new)
    if old_cost and new_cost:
        more = new_cost - old_cost
        percent = more * 100 / old_cost
        if more > 0.0:
            print('Cost has increased by ${0} ({1}%): from ${2} to ${3}'.format(more, percent, old_cost, new_cost))
            sys.exit(0)
        else:
            if more < 0.0:
                down = percent * -1
                print('Cost has decreased by ${0} ({1}%): from ${2} to ${3}'.format(more, down, old_cost, new_cost))
            else:
                print('Cost is the same: ${0}'.format(new_cost))
            print('Everything is awesome!')
            sys.exit(0)
    sys.exit('One or both of the calculated costs is 0!  WTH?: old ({0}) new ({1})'.format(old_cost, new_cost))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ignore_preempted", dest="ignore_preempted", action="store_true", help="ignore preempted tasks")
    parser.add_argument("--only_total", dest="only_total_cost", action="store_true", help="print total cost of the workflow instead of the tsv per task costs")
    either = parser.add_mutually_exclusive_group(required=True)
    either.add_argument('-m', '--metadata', dest='metadata', help='metadata file to calculate cost on')
    either.add_argument('--compare', nargs=2, help='compare old to new cost output')

    args = parser.parse_args()

    if args.metadata:
        with open(args.metadata) as data_file:
            metadata = json.load(data_file)
        calculate_cost(metadata, args.ignore_preempted, args.only_total_cost, True)
        if args.only_total_cost:
            print("Total Cost: " + str(TOTAL_WORKFLOW_COST))
            print("Total run time (hours): " + str(TOTAL_RUN_HOURS))
    else:
        old, new = args.compare
        compare(old, new)


if __name__ == "__main__":
    main()
