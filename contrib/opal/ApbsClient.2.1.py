import sys
import time
import httplib
import string
import os
import subprocess

from AppService_client import \
     AppServiceLocator, getAppMetadataRequest, launchJobRequest, \
     queryStatusRequest, getOutputsRequest, \
     launchJobBlockingRequest, getOutputAsBase64ByNameRequest
from AppService_types import ns0
from ZSI.TC import String

# Help output
if(len(sys.argv)==1 or sys.argv[1]=="-h" or sys.argv[1]=="--help"):
	print open("/home/samir/dev/wustl/opal-py-1.9.3-dev/ApbsClientHelp").read()
	if(len(sys.argv)==2):
		os._exit(99)

# Local run
for arg in sys.argv:
	if(arg=="--local"):
		args=[]
		args.append('apbs')
		for arg2 in sys.argv:
			if(arg2.find("--output-file")!=-1):
				args.append(arg2)
			if(arg2.find("--output-format")!=-1):
				args.append(arg2)
		args.append(sys.argv[-1])
		subprocess.call(args)
		os._exit(99)

# parses custom service location
default=True
for arg in sys.argv:
	if(arg.find("--service-location")!=-1): # change to check beginning of string
		default=False
		if(arg[19:].find("://")!=-1):
			baseURL=arg[19:]
		else:
			baseURL="http://"+arg[19:]
		if(baseURL[-1]!="/"):
			baseURL=baseURL+"/"
		break
	
if(default):
	# Set the protocol to http or https
	proto = "http"
	
	# Host and port for remote services
	#*this string comes from the path to the service
	baseURL = proto + "://ws.nbcr.net:8080/" 


# Retrieve a reference to the AppServicePort
appLocator = AppServiceLocator()
#*this is also from the path to the service
appServicePort = appLocator.getAppServicePort(baseURL + "opal/services/ApbsOpalService") 
	
# Set up remote job launch
req = launchJobRequest()

#*argument parser
# command-line arguments
inFile = sys.argv[-1]
	
# adds any other command line arguments to be passed in req._argList and non-blocking
outFile = False
blocking = True
for arg in sys.argv:
	if(arg.find("--output-file")!=-1): # change to check beginning of string
		req._argList = arg
		outFile=True
	if(arg.find("--output-format")!=-1 and outFile): # change to check beginning of string
		req._argList = req._argList + " " + arg
	if(arg=="--non-blocking"):
		blocking = False

# parses input file
if(inFile.find("/")==-1):
	directory=""
else:
	count=-1
	while inFile[count]!='/':
		count = count-1
	directory = inFile[:count+1]
	inFile = inFile[count+1:]

if(len(sys.argv)==2): # if no other arguments, so it hasn't been created  yet
	req._argList=inFile
else:
	req._argList = req._argList + " " + inFile

#print "argList: \"" + req._argList + "\""

# append all input files in this manner - in this case we have two of them
inputFiles = []
#*this is where apbs.in is read in
inputFiles.append(ns0.InputFileType_Def('inputFile'))
#*this must be the same as req._argList is defined to be
#inputFile0._name = 'apbs.in' 
#print "directory: ", directory
#print "inFile: ", inFile
inputFiles[-1]._name = inFile
tempFile = open(directory+inFile, 'r') 
#tempFileString = tempFile.read()
#inputFiles[-1]._contents = tempFileString
inputFiles[-1]._contents = tempFile.read()
tempFile.close()

# this is where the rest of the files to read in are determined
start = False
tempFile = open(directory+inFile, 'r')
for line in tempFile:
	# remove whitespace
	line=line.strip()
	if(line=="end"):
		break
	if(start and line.find("#")!=0): # eliminates lines with just comments
		# remove comment
		if(line.find("#")!=-1):
			line = line[:line.find("#")]
		# re-remove whitespace (left after comment removal)
		line=line.strip()
		# remove everything except file name
		count = -1
		while line[count]!=' ':
			count = count-1
		fileName=line[count+1:]
		inputFiles.append(ns0.InputFileType_Def('inputFile'))
		inputFiles[-1]._name=fileName
		tempFile2 = open(directory+fileName, "r")
		inputFiles[-1]._contents = tempFile2.read()
		tempFile2.close()
	if(line=="read"):
		start = True

tempFile.close()

#for file in inputFiles:
	#print "fileName: " + file._name
	#if(file._name.find(".in")!=-1):
		#print "contents: " + file._contents


#print "inputFiles: ", inputFiles
# req's inputFile variable is the array of input files created in the lines directly above
req._inputFile = inputFiles

# Launch job, and retrieve job ID
print "Launching remote APBS job"
resp = appServicePort.launchJob(req)
jobID = resp._jobID
print "Received Job ID:", jobID

status = resp._status

if(blocking):
	# Poll for job status
	print "Polling job status"
	while 1:
		# print current status
		print "Status:"
		print "\tCode:", status._code
		print "\tMessage:", status._message
		print "\tOutput Base URL:", status._baseURL

		if (status._code == 8) or (status._code == 4) or (not blocking): # STATUS_DONE || STATUS_FAILED
			break

		# Sleep for 30 seconds
		print "Waiting 30 seconds"
		time.sleep(30)
	    
		# Query job status
		status = appServicePort.queryStatus(queryStatusRequest(jobID))

	# Retrieve job outputs, if execution is successful
	if (status._code==4) or (status._code == 8): # 8 = GramJob.STATUS_DONE
		if(status._code==8):
			print "Retrieving APBS output metadata: "
		else:
			print "The job failed: "
		resp = appServicePort.getOutputs(getOutputsRequest(jobID))

		# Retrieve a listing of all output files
		print "\tStandard Output:", resp._stdOut, "\n", \
			"\tStandard Error:", resp._stdErr
		if (resp._outputFile != None):
			for i in range(0, resp._outputFile.__len__()):
				print "\t" + resp._outputFile[i]._name, ":", resp._outputFile[i]._url
			"\tStandard Error:", resp._stdErr

else:
	print "When the job is complete, the results can be retrieved at: "
	print "\t", status._baseURL

