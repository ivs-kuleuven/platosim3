#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:33:55 2022

@author: juano

adapted from
https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access
"""

#ASYNCHRONOUS REQUEST

#Python 2
#import httplib
#import urllib

#Python 3
import http.client as httplib
import urllib.parse as urllib
import time
from xml.dom.minidom import parseString

host = "gea.esac.esa.int"
port = 443
pathinfo = "/tap-server/tap/async"

import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io.votable import parse

import pandas as pd

#%% Configuration variables
lopn1 = SkyCoord( ra = 277.18023*u.deg, dec =  52.85952*u.deg, frame = 'icrs')
lops1 = SkyCoord( ra =  93.49134*u.deg, dec = -42.93544*u.deg, frame = 'icrs')
plato_fov = 18.*u.deg

catalogue = 'gaiadr3.gaia_source'
maglim    = 10

outputFileName = "fgs_example.vot"


#%% Create job
#-------------------------------------
#Create job

coord  = lops1.copy()
radius = plato_fov

params = urllib.urlencode({\
	"REQUEST"        : "doQuery",       \
	"LANG"           : "ADQL",          \
	"FORMAT"         : "votable_plain", \
	"PHASE"          : "RUN",           \
	"JOBNAME"        : "PLATO FGS catalog", \
	"JOBDESCRIPTION" : "Masterarbeit S. Bowling (contact juan.cabrera@dlr.de)", \
	"QUERY"          : f"SELECT DISTANCE(POINT({coord.ra.deg},{coord.dec.deg}),POINT(ra,dec)) AS dist, designation, ra, dec, phot_g_mean_mag FROM {catalogue} AS cat WHERE 1=CONTAINS(POINT({coord.ra.deg},{coord.dec.deg}),CIRCLE(cat.ra,cat.dec,{radius.value})) AND cat.phot_g_mean_mag < {maglim} ORDER BY dist ASC"
	})

    
headers = {\
	"Content-type": "application/x-www-form-urlencoded", \
	"Accept":       "text/plain" \
	}

connection = httplib.HTTPSConnection(host, port)
connection.request("POST",pathinfo,params,headers)

#Status
response = connection.getresponse()
print ("Status: " +str(response.status), "Reason: " + str(response.reason))

#Server job location (URL)
location = response.getheader("location")
print ("Location: " + location)

#Jobid
jobid = location[location.rfind('/')+1:]
print ("Job id: " + jobid)

connection.close()

#-------------------------------------
#Check job status, wait until finished

while True:
	connection = httplib.HTTPSConnection(host, port)
	connection.request("GET",pathinfo+"/"+jobid)
	response = connection.getresponse()
	data = response.read()
	#XML response: parse it to obtain the current status
	#(you may use pathinfo/jobid/phase entry point to avoid XML parsing)
	dom = parseString(data)
	phaseElement = dom.getElementsByTagName('uws:phase')[0]
	phaseValueElement = phaseElement.firstChild
	phase = phaseValueElement.toxml()
	print ("Status: " + phase)
	#Check finished
	if phase == 'COMPLETED': break
	#wait and repeat
	time.sleep(0.2)


connection.close()

#-------------------------------------
#Get results
connection = httplib.HTTPSConnection(host, port)
connection.request("GET",pathinfo+"/"+jobid+"/results/result")
response = connection.getresponse()
data = response.read().decode('iso-8859-1')
#print(type(data))
#print(data)
#outputFileName = "example_votable_output.vot"
outputFile = open(outputFileName, "w")
outputFile.write(data)
outputFile.close()
connection.close()
print ("Data saved in: " + outputFileName)	

#%%

# https://gist.github.com/icshih/52ca49eb218a2d5b660ee4a653301b2b
def votable_to_pandas(votable_file):
    votable = parse(votable_file)
    table = votable.get_first_table().to_table(use_names_over_ids=True)
    return table.to_pandas()

tb = votable_to_pandas( outputFileName)
tb.to_csv( outputFileName.replace( '.vot', '.csv'))
