import sys
import time
import httplib

from AppService_client import \
     AppServiceLocator, getAppMetadataRequest, launchJobRequest, \
     queryStatusRequest, getOutputsRequest, \
     launchJobBlockingRequest, getOutputAsBase64ByNameRequest
from AppService_types import ns0
from ZSI.TC import String

# Set the protocol to http or https
proto = "http"

# Host and port for remote services
baseURL = proto + "://ws.nbcr.net:8080/"

# Retrieve a reference to the AppServicePort
appLocator = AppServiceLocator()
appServicePort = appLocator.getAppServicePort(
    baseURL + "opal/services/ApbsOpalService")
	
# Set up remote job launch
req = launchJobRequest()

# command-line arguments
req._argList = "apbs.in"

# append all input files in this manner - in this case we have two of them
inputFiles = []
inputFile0 = ns0.InputFileType_Def('inputFile')
inputFile0._name = 'apbs.in'
sampleFile0 = open("etc/apbs.in", "r")
sampleFileString0 = sampleFile0.read()
sampleFile0.close()
inputFile0._contents = sampleFileString0
inputFiles.append(inputFile0)

inputFile1 = ns0.InputFileType_Def('inputFile')
inputFile1._name = 'ion.xml'
sampleFile1 = open("etc/ion.xml", "r")
sampleFileString1 = sampleFile1.read()
sampleFile1.close()
inputFile1._contents = sampleFileString1
inputFiles.append(inputFile1)

req._inputFile = inputFiles

# Launch job, and retrieve job ID
print "Launching remote Pdb2pqr job"
resp = appServicePort.launchJob(req)
jobID = resp._jobID
print "Received Job ID:", jobID

# Poll for job status
status = resp._status
print "Polling job status"
while 1:
    # print current status
    print "Status:"
    print "\tCode:", status._code
    print "\tMessage:", status._message
    print "\tOutput Base URL:", status._baseURL

    if (status._code == 8) or (status._code == 4): # STATUS_DONE || STATUS_FAILED
        break

    # Sleep for 30 seconds
    print "Waiting 30 seconds"
    time.sleep(30)
    
    # Query job status
    status = appServicePort.queryStatus(queryStatusRequest(jobID))

# Retrieve job outputs, if execution is successful
if status._code == 8: # 8 = GramJob.STATUS_DONE
    print "Retrieving Pdb2pqr output metadata: "
    resp = appServicePort.getOutputs(getOutputsRequest(jobID))

    # Retrieve a listing of all output files
    print "\tStandard Output:", resp._stdOut, "\n", \
	  "\tStandard Error:", resp._stdErr
    if (resp._outputFile != None):
	for i in range(0, resp._outputFile.__len__()):
	    print "\t" + resp._outputFile[i]._name, ":", resp._outputFile[i]._url
