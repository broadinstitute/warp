from firecloud import fiss
import firecloud.api as fapi
import os
import io
import pandas as pd
import argparse
import sys

billing_project='hca-dev-staging-billing'
template_workspace_name="DCP2_Optimus_template_FK"


def create_newworkspace(billing_project, template_workspace_name, new_workspace_name):

   ent_types = fiss.fapi.list_entity_types(billing_project, template_workspace_name)

   res = fiss.fapi.clone_workspace(billing_project, 
                                   template_workspace_name, 
                                   billing_project, 
                                   new_workspace_name
                                  )

   print(res.text)


def main():
   global billing_project
   global template_workspace_name

   parser = argparse.ArgumentParser(prog = "python " + sys.argv[0], add_help = False)
   subparser = parser.add_subparsers(dest="cmd")

   delete_workspace = subparser.add_parser('delete_workspace', help = 'delete workspace')
   delete_workspace.add_argument('--name', help = "name of the workspace")

   clone_workspace = subparser.add_parser('clone_workspace',   help = 'clone from existing workspace')
   clone_workspace.add_argument('--source-work-space', dest = 'src_work_space',  help = "name of source workspace")
   clone_workspace.add_argument('--destination-work-space', dest = 'dest_work_space', help = "name of destination workspace")

   # show help when no arguments supplied
   if len(sys.argv) == 1:
       parser.print_help()
       sys.exit(0)

   args = parser.parse_args()

   #new_workspace_name = "DCP2_Optimus_template_KMK_v1"
   if args.cmd == 'delete_workspace':
      print("Delete existing workspace ", args.name)
      delete_status =  fapi.delete_workspace(billing_project, args.name)

   elif args.cmd == 'clone_workspace':
      print("Cloning a new workspace from template", args.dest_work_space)
      status = create_newworkspace(billing_project, args.src_work_space, args.dest_work_space)

if __name__ == "__main__":
    main()

