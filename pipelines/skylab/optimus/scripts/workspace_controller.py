from firecloud import fiss
import firecloud.api as fapi
import os
import io
import pandas as pd
import argparse
import sys
import parse_terra

billing_project= os.environ['WORKSPACE_NAMESPACE']
template_workspace_name="DCP2_Optimus_template_FK"


def create_newworkspace(billing_project, template_workspace_name, new_workspace_name):

   ent_types = fiss.fapi.list_entity_types(billing_project, template_workspace_name)

   res = fiss.fapi.clone_workspace(billing_project, 
                                   template_workspace_name, 
                                   billing_project, 
                                   new_workspace_name
                                  )

   print(res.text)

def upload_tables(input_file, billing_project, workspace_name):
      with open(input_file) as tsvf:
         headerline = tsvf.readline().strip()
         entity_data = [l.rstrip('\n') for l in tsvf]

      model = 'firecloud'
      fiss._batch_load(billing_project, workspace_name, headerline, entity_data, 500, model)




def main():
   global billing_project
   global template_workspace_name

   parser = argparse.ArgumentParser(prog = "python " + sys.argv[0], add_help = False)
   subparser = parser.add_subparsers(dest="cmd")

   delete_workspace = subparser.add_parser('delete_workspace', help = 'delete workspace')
   delete_workspace.add_argument('--workspace-name', dest = "workspace_name",  help = "name of the workspace")

   clone_workspace = subparser.add_parser('clone_workspace',   help = 'clone from existing workspace')
   clone_workspace.add_argument('--source-work-space', dest = 'src_work_space',  help = "name of source workspace")
   clone_workspace.add_argument('--destination-work-space', dest = 'dest_work_space', help = "name of destination workspace")

   get_data_info = subparser.add_parser('get_data_info',   help = 'get participant.tsv')
   get_data_info.add_argument('--workspace-name', dest="workspace_name",  help = "name of the workspace")
   get_data_info.add_argument('--participant-table-name', dest="participant_table_name",  help = "name of sample table")
   get_data_info.add_argument('--output-name', dest="output_table_name", required = False, 
                               default = "participant.tsv",  help = "name of output tsv"
                             )

   create_participant_lane = subparser.add_parser('create_participant_lane',
                                                  help = 'create participant_lane/lane_set_id tables')
   create_participant_lane.add_argument('--input-name', dest="input_participant_table_name", 
                                        required = False, 
                                        default = "participant.tsv",  help = "input participant table  name")

   create_participant_lane.add_argument('--output-prefix', dest="output_prefix", required = False, 
                               help = "name of output prefix for the lanes")


   upload_participant_lane = subparser.add_parser('upload_participant', 
                                                  help = 'uploads the participant_lane_set, _lane_membership and _lane_entity files')
   upload_participant_lane.add_argument('--workspace-name', dest="workspace_name",  help = "name of the workspace")
   upload_participant_lane.add_argument('--input-prefix', dest="input_prefix",  help = "name of the input prefix")

   # show help when no arguments supplied
   if len(sys.argv) == 1:
       parser.print_help()
       sys.exit(0)

   args = parser.parse_args()

   #new_workspace_name = "DCP2_Optimus_template_KMK_v1"
   if args.cmd == 'delete_workspace':
      print("Delete existing workspace ", args.workspace_name)
      delete_status =  fapi.delete_workspace(billing_project, args.workspace_name)

   elif args.cmd == 'clone_workspace':
      print("Cloning a new workspace from template", args.dest_work_space)
      status = create_newworkspace(billing_project, args.src_work_space, args.dest_work_space)

   elif args.cmd == 'get_data_info':
      print("Get information from workspace", args.workspace_name)
      r = fapi.get_entities_tsv(billing_project, args.workspace_name, args.participant_table_name)
      with open(args.output_table_name, 'w') as fout: 
          fout.write(r.content.decode())

   elif args.cmd == 'create_participant_lane':
       parse_terra.create_output_files(args.input_participant_table_name, args.output_prefix)

   elif args.cmd == 'upload_participant':
      #upload_tables(args.input_prefix + ".tsv", billing_project, args.workspace_name)
      upload_tables(args.input_prefix + "_lane.tsv", billing_project, args.workspace_name)
      #upload_tables(args.input_prefix + "_lane_membership.tsv", billing_project, args.workspace_name)
      #upload_tables(args.input_prefix + "_lane_entity.tsv", billing_project, args.workspace_name)


if __name__ == "__main__":
    main()

