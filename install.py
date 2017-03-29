#!/usr/bin/env python

from __future__ import print_function


__author__ = "Joris van Steenbrugge"

import os
import subprocess as sp
from sys import stderr
from glob import glob




def parseConfig(path):   
    path = path.replace("/","\/")
    sp.call("sed -i 's/ac_default_prefix=\/usr\/local/ac_default_prefix={}/g' configure".format(
                    path), 
                shell = True)

def checkOutfolder(path):
    content = glob(path+"/*")
    for c in content:
        if "bin" in c:
            sp.call("mv {} {}".format(path+"/bin/*", path), 
                        shell = True)
            sp.call("rm -r {}".format(path+"/bin"),
                        shell = True)

def downloadtmp(url, name):
    os.chdir("/tmp")

    print("Downloading source code to /tmp/", file=stderr)

    sp.check_output("wget {} -O prog.tar.gz".format(url),
            stderr=sp.STDOUT, shell=True)
    
    sp.call("mkdir progtar && tar xf prog.tar.gz -C progtar --strip-components 1",
                shell=True)

    os.chdir("progtar")

    try:

        if os.path.exists("configure"):
            print("Configure file was found.. changing the install directory", file=stderr)
            parseConfig(wd+name)
            print("Running configure", file=stderr)
            sp.check_output("./configure", 
                   shell=True, stderr=sp.STDOUT)
            print("Compiling the source", file=stderr)
            sp.check_output("make  && make install ", 
                    shell=True, stderr=sp.STDOUT)

        elif os.path.exists("Makefile") or os.path.exists("makefile"):
            pass #not implemented
        else:
            sp.call("mv * {}".format(wd+name),
                    shell=True)
    except sp.CalledProcessError:
        pass

    checkOutfolder(wd+name)
    sp.call("rm /tmp/progtar -r", shell=True)    

def setup(tools):
    wd = os.getcwd()+"/tools/"

    print("Creating tools path: {}".format(wd), file=stderr)
    try:
        os.mkdir(wd)
    except OSError:
        print("Path already exists", file=stderr)
    

    print("Tools will be installed in:", file=stderr)
    for tool in tools.keys():
        try:
            path = wd+tool
            os.mkdir(path)
            print("    "+path, file=stderr)
        except OSError:
            print("    {} seems to be already installed".format(tool))
            del tools[tool]
    return wd, tools


if __name__ == "__main__":
    tools = {
            "idba" : "https://github.com/loneknightpy/idba/releases/download/1.1.3/idba-1.1.3.tar.gz",
            "MetaSPAdes": "http://cab.spbu.ru/files/release3.10.1/SPAdes-3.10.1-Linux.tar.gz",
            "megaHit": "https://github.com/voutcn/megahit/releases/download/v1.1.1/megahit_v1.1.1_LINUX_CPUONLY_x86_64-bin.tar.gz"
            }

    wd, tools = setup(tools)

    for tool, url in tools.items():
        downloadtmp(url, tool)

