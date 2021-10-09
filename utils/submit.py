
#===================================
#        Chengwen Liu              #
#      liuchw2010@gmail.com        #
#   University of Texas at Austin  #
#===================================

import os
import sys
import argparse
import subprocess

def subone(node, command):
  subcmd = f'ssh {node} "cd {path}; export CUDA_DEVICE_ORDER=PCI_BUS_ID; export CUDA_VISIBLE_DEVICES=0; {command}"'
  os.system(subcmd)
  return 

def main():
  parser = argparse.ArgumentParser()
  parser.add_argument('-n', dest = 'node', help = "node to submit", required = True)  
  parser.add_argument('-p', dest = 'path', help = "path to submit", default = ' ')  
  parser.add_argument('-c', dest = 'command', nargs='+', help = "command strings", default = [], required=True)  
  args = vars(parser.parse_args())
  global path
  node = args["node"]
  path = args["path"]
  command = args["command"]
  
  if path == ' ':
    path = os.getcwd()
  else:
    path = os.path.join(os.getcwd(), path)
  command = '  '.join(command)
  rtx30s = ['node36', 'node102', 'bme-mercury', 'mercury', 'bme-rna', 'rna']
  if node in rtx30s:
    if 'cuda10.2' in command:
      command = command.replace('cuda10.2', 'cuda11.2')
  subone(node, command)
  return

if __name__ == "__main__":
  main()
