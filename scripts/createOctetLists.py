#!/usr/bin/python

import os

if 0:
	with open("BetaRunsInOctet_2011-2012.txt", "wt") as out:
		with open("OctetList_20112012.txt") as fp:
			for line in fp:
				if (line[:6]=="Octet:"):
					words=line.split()
					for idx, w in enumerate(words):
						if w=='A2' or w=='A5' or w=='A7' or w=='A10' or w=='B2' or w=='B5' or w=='B7' or w=='B10': 
							out.write(words[idx+2] + ' ')
					out.write('\n')

if 1:
	with open("BetaRunsInOctet_2012-2013.txt", "wt") as out:
		with open("OctetList_20122013.txt") as fp:
			for line in fp:
				if (line[:6]=="Octet:"):
					words=line.split()
					for idx, w in enumerate(words):
						if w=='A2' or w=='A5' or w=='A7' or w=='A10' or w=='B2' or w=='B5' or w=='B7' or w=='B10': 
							out.write(words[idx+2] + ' ')
					out.write('\n')
				


octet = 1
				
with open(os.getenv("OCTET_LIST")+"OctetList_20112012.txt") as fp:
	for line in fp:
		if (line[:6]=="Octet:"):
			with open(os.getenv("OCTET_LIST")+"octet_list_%i.dat"%octet, "wt") as out:
				words=line.split()
				for idx, w in enumerate(words):
					if w[:1]=='A' or w[:1]=='B':
						out.write(words[idx] + ' ' + words[idx+2] + '\n')
				octet+=1
