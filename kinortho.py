#-----------------------------------------------------------------------------
#----------------------------NEEDS TO RUN WITH LINUX--------------------------
#-----------------------------------------------------------------------------
import sys
import Bio          #------------ importing stuff----------------------------
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import os
import math
#--------------------------------NEED BIOPYTHON TO RUN!-----------------------
#--------------------------------NEED BLAST+ TO RUN !-------------------------
#--------------------------------NEED NUMPY to RUN !!!!!!_--------------------


#-----------------------------------Helper varibal/input stuff----------------
infolder = sys.argv[1]
x=2
fullseq=""
domseq=""
out=""
evalue=""
thread=""
inflat=""
minev=""

#woop="1234567"
#woop=woop[0:8]
#print (woop)

#-----------------just taking the input from bash and making------------------
#-----------------it be able to take any comand in any order------------------

for i in range(int(len(sys.argv)/2)-1):#loops in pairs through the bash command
    if len(sys.argv)>x:                #       if the commands are even present
        if sys.argv[x]=="-f":
            fullseq = sys.argv[x+1]
        if sys.argv[x]=="-d":
            domseq = sys.argv[x+1]
        if sys.argv[x]=="-o":
            out = sys.argv[x+1]
        if sys.argv[x]=="-E":
            evalue = sys.argv[x+1]  #         all these conditionals check what
        if sys.argv[x]=="-t":       #              are present in the bash line
            thread = sys.argv[x+1]
        if sys.argv[x]=="-i":
            inflat = sys.argv[x+1]
        if sys.argv[x]=="-e":
            minev = sys.argv[x+1]
    x=x+2
if minev == "":
    minev = float(1e-200)
if evalue == "":
    evalue=float(1e-5)
if inflat== "":
    inflat=1.5
#firstcommand= "makeblastdb -in ref.fasta -dbtype prot -out ref"

#-------------------------------defining some methods/funtions-----------------
def build_blast_db(input,dbname):
    firstcommand= "makeblastdb -in "
    firstcommand= firstcommand + input
    firstcommand= firstcommand +" -dbtype prot -out "
    firstcommand= firstcommand + dbname
    os.system(firstcommand)


def homology_search(inf,db,outf):
    command="blastp -outfmt \"6 qseqid sseqid sstart send evalue bitscore positive\""
    command= command + " -query " + inf
    command= command + " -db " + db
    if(thread!=""):
        command= command + " -num_threads " + thread

    if(evalue!=""):
        command= command + " -evalue " + str(evalue)

    # add condtionals
    command= command + " > " + outf
    os.system(command)


class Edgy:
    def __init__(self,location,weight):
        self.w = weight
        self.loc = location

    def getW(self):
        return self.w

    def getL(self):
        return self.loc
#--------------------------------------------------------------------------

index=0  # helper varibal
rec=[]
tax=[]
amountofspec=0

records = list(SeqIO.parse(fullseq, "fasta"))
fullrec=[]
for record in records:
    r=record.id
    #r=r+'#' + infile
    record.id = r
    fullrec.append(r)
    index= index +1
#print(fullrec)

#-------------------------This is for retreaving files------------------------
for infile in os.listdir(infolder):
    amountofspec=amountofspec+1
    goodfile = infolder + infile
    records = list(SeqIO.parse(goodfile, "fasta"))      # reads the input files
    for record in records:
        r=record.id
        r=r+'#' + infile
        record.id = r
        rec.append(record)
        index= index +1
print(index)
#---------------------------This is for writing a file------------------------
ref=open("ref.fasta","w+")
SeqIO.write(rec,"ref.fasta","fasta")
ref.close()





#woop=Edgy("asdklfal",3)
#print(woop.getW())
#print(woop.getL())


#'''

#-------------First step in the program to build a database--------------------
build_blast_db("ref.fasta","firstDB")
#---------------Blastp the full seq agasit the database------------------------
#'''
#'''

homology_search(fullseq,"firstDB","txtHelper.txt")

print("first homology search done")


helperfile=open("txtHelper.txt","r")
lines = helperfile.readlines()

#-----------------------read the results form the Blastp-----------------------
helperfile=open("txtHelper.txt","r")
lines = helperfile.readlines()

sseqids=[]
tax={}
fulltax={}



curfull="woop"

for line in lines:
    b=0
    b2=0
    fcur=""
    fkey=""
    ftax=""
    for char in line:
        if char=='\t':
            fkey=fcur
        if b == 0:
            if char != '\t':
                fcur=fcur+char
        else:
            if char=='#':
                b2=1

            if b2==1:
                if char=='\t':
                    break
                if char != '\t':
                    ftax=ftax+char
        if char=='\t':
            b=1
    if curfull!=fkey:
        curfull=fkey
        fulltax[fkey]=ftax


#print(fulltax)

index=0

for line in lines:  #                for each line in the result in the blastp
    l=list(line)
    index=0
    tempstring=""
    helperbool=0
    b=0
    for chara in l:  #                    for each chararicter in the line
        if chara == '\t':
            for i in range(len(l)):#        record teh line until the tab char
                if l[index+i+1] != '\t':
                    tempstring= tempstring + l[index+i+1]
                else:
                    for sseqid in sseqids: #                   for each sseqid
                        if str(sseqid) == tempstring:
                            helperbool=1
                    if helperbool ==0:
                        sseqids.append(tempstring)#add ids to the array of ids
                        #print(tempstring)
                    break
            break
        index=index+1



#-----------------------geting the files ready for all vs all blast------------

newrec=[]
for record in rec:


    for sseqid in sseqids:
        if record.id==sseqid:
            newrec.append(record)

print(len(rec))
print(len(sseqids))
print(len(newrec))

#----------------read and writes the file for the second blastP---------------
newref=open("homology.fasta","w+")
SeqIO.write(newrec,"homology.fasta","fasta")
newref.close()

build_blast_db("homology.fasta","homology")

#----------------blastp of all vs all----------------------------------------
print("allvsall start")

homology_search("homology.fasta","homology","all_vs_all.blastp")


print("allvsall finish")


#----------------------------------------------------------------------------
#----------reads and retrives the data from the output of the----------------
#------------all vs all blastP and puts into dictionaries--------------------
#----------------------------------------------------------------------------

helperfile=open("all_vs_all.blastp","r")
lines = helperfile.readlines()

keys=[]
halfkeys=[]
indhalfkeys=[]
keysE={}
keysB={}
keysP={}
minE={}
maxB={}
maxP={}
tax={}
seckeys={}
firkeys={}
samebool=[]
taxorder={}
#sstart={}
#send={}
taxorderhelper=0

#aveP={}
#-------------go through the files and find all the info of the -------------
#--------------------AllvsAll Blastp outputfile------------------------------

#          things you should know about about blastp files
#          there differnt info is sesperated by a tab (\t)
#       in the code when you see tab(\t) its the end thing inbetween info


for line in lines:# for each line
    l=list(line)
    #print(l)
    #helper varibles
    index=0
    fin =0
    f=""
    s=""
    firststring=""
    tempstring=""
    helperbool=0
    taxbool=0
    taxa=""

    for chara in l: #                      for each character in each line
        if chara !='\t':#     if not the end of the line add to temp sting
            tempstring= tempstring +chara
            if helperbool == 0:#          find first string and add to var
                firststring= firststring+ chara
                f=f+chara
            else:
                s=s+chara
            #print(tempstring)
        if chara == '\t':#   if end of string in file add taxonomy to dic
            taxbool=0
            tax[tempstring]=taxa
            if not taxa in taxorder.keys():#       checks for no reapeats
                taxorder[taxa]=taxorderhelper
                taxorderhelper=taxorderhelper+1
            if helperbool==1:#               if the second id was passed
                if f == s: #      if first id and sencond id is the same
                    samebool.append(1)#to be honest i dont even use this
                else: #         but to scared to scared to get rid of it
                    samebool.append(0)
                index=index+1       #index helps counts the where we are
                keys.append(tempstring) #    add the full key to the dic
                charbool=0              #       We are done with id part
                secstring=""
                for char in tempstring: #   for the charartors in stirng
                    if charbool==1: #    if on second id place it in dic
                        secstring=secstring+char
                    if char=='\t':
                        charbool=1
                firkeys[tempstring]=firststring
                seckeys[tempstring]=secstring

                if not firststring in indhalfkeys:#checking for repeats
                    indhalfkeys.append(firststring)
#                      make a array of half keys to help later funtions
#                                       Both unquie and non unquie list
                halfkeys.append(firststring)
                #if not tempstring in aveB.keys():
                #    aveB[tempstring]=[]
                #    aveP[tempstring]=[]

        #------ after the ID and into the other data------
                #ss=""
                #if secstring not in sstart.keys():
                #    sstart[secstring]=[]
                for i in range(len(l)):# the start info and redoring it
                    if l[index]!='\t':
                #        ss=ss+l[index]
                        index=index+1
                    if l[index]=='\t':
                        index=index+1
                        break
                #if ss not in sstart[secstring]:
                #sstart[secstring].append(ss)
                #se=""
                #if secstring not in send.keys():
                #    send[secstring]=[]
                for i in range(len(l)):#  the end info and recording it
                    if l[index]!='\t':
                #        se=se+l[index]
                        index=index+1
                    if l[index]=='\t':
                        index=index+1
                        break
                #if se not in send[secstring]:
                #send[secstring].append(se)
                e=""#                      the E-value and recording it
                for i in range(len(l)):
                    if l[index]!='\t':
                        e=e+l[index]
                        index=index+1
                    if l[index]=='\t':
                        keysE[tempstring]=e
                        if not firststring in minE.keys():
                            minE[firststring]=e
                        else:
                            if float(e) < float(minE[firststring]):
                                minE[firststring]=e
                        index=index+1
                        break
                b=""#                   the bit-score and recording it
                for i in range(len(l)):
                    if l[index]!='\t':
                        b=b+l[index]
                        index=index+1
                    if l[index]=='\t':
                        keysB[tempstring]=b
                        if not firststring in maxB.keys():
                            maxB[firststring]=b
                        else:
                            if float(b) > float(maxB[firststring]):
                                maxB[firststring]=b

                        index=index+1
                        break
                p=""#                    The postive and recording it
                for i in range(len(l)):#   for the length of the line
                    if l[index]!='\t':#            tab inbetween info
                        p=p+l[index]
                        index=index+1 #             add more to index
                    if index==len(l):#                 if end of line
                        keysP[tempstring]=p#
                        if not firststring in maxP.keys():
                            maxP[firststring]=p#       add pos to dic
                        else:
                            if float(p) > float(maxP[firststring]):
                                maxP[firststring]=p
                        fin=1
                        # helper bool that tells us that we are done
                        index=index+1 #          gotta add the index
                        break
            tempstring=tempstring+'\t'
            helperbool=1
            taxa=taxa+'#'

        if fin == 1:#                                  if done break
            break

        if taxbool==1:     #                   stuff for the tax dic
            taxa=taxa+chara
        if chara=='#':
            taxbool =1
        index=index+1

helperfile.close()


#print(len(keys))
#for key in keys:
#    print(key)

print("allvall scanner all done")
#----------------------creats helpers comp sturctures-----------------------

tax2={}#                               helps dics for the taxonomy
tax1={}
taxorder={}
taxorderhelper=0

for key in keys: #                                for the all keys
    t1=""#                                               first tax
    t2=""#                                              second tax
    taxa=tax[key]#      gets the tax dic and slits it into two dic
    helperbool=0
    for char in taxa:
        if char=='#':
            helperbool=1
            continue
        if helperbool == 0:
            t1=t1+char
        if helperbool ==1:
            t2=t2+char
    if not t2 in taxorder.keys():
        taxorder[t2]=taxorderhelper
        taxorderhelper=taxorderhelper+1
    tax1[key]=t1
    tax2[key]=t2

#-----------------------finding paralogs------------------------------------
para={}
i=0
cur="woop"
helperbool=0
helperbool2=0
for key in keys:

    if cur != halfkeys[i]:
        cur=halfkeys[i]
        helperbool=0
        para[halfkeys[i]]=[]
        helperbool2=0


    if tax2[key]!=tax1[key]:
        helperbool=1
    if helperbool == 0:
        if helperbool2==1:
            #print(key)
            #print(keysE[key])
            if float(keysE[key]) <= evalue:
                para[halfkeys[i]].append(key)
            #    print(key)
    helperbool2=1
    i=i+1

#----------------------------finding orthologs---------------------------

ortho={}
i=0
cur="woop"
helperbool=0
helperbool2=0
boollist=[]
for x in range(amountofspec):
    boollist.append(0)
for key in keys:

    if cur != halfkeys[i]:
        cur=halfkeys[i]
        #helperbool=0
        ortho[halfkeys[i]]=[]
        for x in range(amountofspec):
            ortho[halfkeys[i]].append("NA")
            boollist[x]=0
        #print(key)
        #print(tax2[key])
        #print(taxorder[tax2[key]])
        boollist[taxorder[tax2[key]]]=1


    if boollist[taxorder[tax2[key]]]==0:
        if float(keysE[key]) <= evalue:
            #print(key)
            ortho[halfkeys[i]][taxorder[tax2[key]]]=key

        #    print(tax2[key])
        #    print(taxorder[tax2[key]])
            boollist[taxorder[tax2[key]]]=1

    i=i+1


#--------------------------finding co-ortholongs ---------------------------

coortho={}

for key in indhalfkeys:
    coortho[key]=[]
    for okey in ortho[key]:
        if okey!="NA":
            if seckeys[okey] in para.keys():
                for pkey in para[seckeys[okey]]:
                    coortho[key].append(pkey)
                    #print(pkey)
#-------------------spitting out the output---------------------------------


loglist=[0,0]
indexlist=[0,0]
for key in indhalfkeys:
    for orth in ortho[key]:
        if orth!="NA":
            #print(orth)
            #print(keysE[orth])
            if float(keysE[orth]) != 0.0:
                e=-(math.log10(float(keysE[orth])))

                loglist[0]=loglist[0]+e
            else:
                loglist[0]=loglist[0]+200
                e=200
            indexlist[0]=indexlist[0]+1
            #print(e)
    for corth in coortho[key]:

        #print(corth)
        if float(keysE[corth]) != 0.0:
            e=-(math.log10(float(keysE[corth])))
        #print(e)
            loglist[0]=loglist[0]+e
        else:
            loglist[0]=loglist[0]+200
            e=200
        indexlist[0]=indexlist[0]+1
    for par in para[key]:

        #print(corth)
        if float(keysE[par]) != 0.0:
            e=-(math.log10(float(keysE[par])))
        #print(e)
            loglist[1]=loglist[1]+e
        else:
            loglist[1]=loglist[1]+200
            e=200
        indexlist[1]=indexlist[1]+1
#print(loglist[0])
#print(indexlist[0])
normortho=loglist[0]/indexlist[0]
normpara=loglist[1]/indexlist[1]



output=open("graph.abc", "w+")

network={}
for key in indhalfkeys:
    vert=[]
    for orth in ortho[key]:
        if orth!="NA":
            if float(keysE[orth]) != 0.0:
                e=-(math.log10(float(keysE[orth])))
                e=e/normortho
            else:
                e=200/normortho

            woop=str(orth+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(orth,e)
            vert.append(edge)
            #print(edge.getL())
            #print(edge.getW())
    for corth in coortho[key]:
        if float(keysE[corth]) != 0.0:
            e=-(math.log10(float(keysE[corth])))
            e=e/normortho
        else:
            e=200/normortho
        woop=str(corth+'\t'+str(e))
        output.write(woop)
        output.write("\n")
        edge=Edgy(corth,e)
        vert.append(edge)
        #print(edge.getL())
        #print(edge.getW())
    for par in para[key]:
        if float(keysE[par]) != 0.0:
            e=-(math.log10(float(keysE[par])))
            e=e/normpara
        else:
            e=200/normpara
        woop=str(par+'\t'+str(e))
        output.write(woop)
        output.write("\n")
        edge=Edgy(par,e)
        vert.append(edge)
        #print(edge.getL())
        #print(edge.getW())

    network[key]=vert

output.close()

os.system("mcxload -abc graph.abc -o OUTPUT.mci -write-tab OUTPUT.tab")

com = "mcl OUTPUT.mci -I "
com = com + str(inflat)
com = com + " -o OUTPUT.out"

os.system(com)
#--------------------

#---------------------
tab2prot={}
keykeys=[]

helperfile=open("OUTPUT.tab","r")
lines = helperfile.readlines()
for line in lines:
    prot=""
    tab=""
    helperbool=0
    for char in line:
        if char == '\t':
            helperbool=1
            continue
        if helperbool==0:
            tab=tab+char
        else:
            if char=='\n':
                continue
            prot=prot+char
    #print(tab)
    #print(prot)
    tab2prot[tab]=prot
    keykeys.append(prot)


prot2cluster={}




helperfile.close()

helperfile=open("OUTPUT.out","r")
lines = helperfile.readlines()
x=0
prot=""
cluster=""
prot2clus={}
for line in lines:
    if x<7:
        x=x+1
        continue
    prot=""
    helperbool=0
    hb=0
    if list(line)[0]!=' ':
        #prot=""
        cluster=""
        #helperbool=0
        #hb=0
        if list(line)[0] ==')':
            break
        for char in line:
            if char==" " or char=='\n':
                if helperbool==1:
                    if hb==1:
                        if prot != "$":
                            if tab2prot[prot] not in prot2clus.keys():
                                prot2clus[tab2prot[prot]]=[]
                            prot2clus[tab2prot[prot]].append(cluster)
                            #print(tab2prot[prot])
                            #print(cluster)
                prot=""
                helperbool=1
                continue
            if helperbool==0:
                cluster=cluster+char
            else:
                hb=1
                prot=prot+char
    else:
        x=0
        prot=""
        for char in line:

            if x<7:
                x=x+1
                continue

            if char==" " or char=='\n':
                if prot != "$":
                    if tab2prot[prot] not in prot2clus.keys():
                        prot2clus[tab2prot[prot]]=[]
                    prot2clus[tab2prot[prot]].append(cluster)
                    #print(tab2prot[prot])
                    #print(cluster)
                    prot=""
            else:
                prot=prot+char

goodpara=[]
goodorth=[]
goodco=[]
i1=0
i2=0
for key in indhalfkeys:

    for pkey in para[key]:
        if firkeys[pkey]!=seckeys[pkey]:
            i1=i1+1
            if firkeys[pkey] in prot2clus.keys() and seckeys[pkey] in prot2clus.keys():
                for x in prot2clus[firkeys[pkey]]:
                    for y in prot2clus[seckeys[pkey]]:
                        if x == y:
#                        print(x)
#                        print(y)
#                        print(pkey)
                            goodpara.append(pkey)
                            i2=i2+1
    for ckey in coortho[key]:
        if firkeys[ckey]!=seckeys[ckey]:
            i1=i1+1
            if firkeys[ckey] in prot2clus.keys() and seckeys[ckey] in prot2clus.keys():
                for x in prot2clus[firkeys[ckey]]:
                    for y in prot2clus[seckeys[ckey]]:
                        if x == y:
                            goodco.append(ckey)
                            i2=i2+1


    for okey in ortho[key]:
        if okey != "NA" :
            i1=i1+1
            if firkeys[okey] in prot2clus.keys() and seckeys[okey] in prot2clus.keys():
                for x in prot2clus[firkeys[okey]]:
                    for y in prot2clus[seckeys[okey]]:
                        if x == y:
                            goodorth.append(okey)
                            i2=i2+1



print(i1)
print(i2)

tempstr=""
templist=[]

for key in fullrec:
    tempstr=""
    tempstr=key+fulltax[key]
    templist.append(tempstr)
fullrec=templist
#print(fullrec)

goodclus=[]

for key in fullrec:
    if key in prot2clus.keys():
        for clus in prot2clus[key]:
            if clus not in goodclus:

                goodclus.append(clus)

#print(goodclus)
forth=[]
fpara=[]
fco=[]
f1=0
f2=0
f3=0
for orth in goodorth:
    for cluster in goodclus:
        for y in prot2clus[seckeys[orth]]:
            if y==cluster:
                forth.append(orth)
                f1=f1+1
for para in goodpara:
    for cluster in goodclus:
        for y in prot2clus[seckeys[para]]:
            if y==cluster:
                fpara.append(para)
                f2=f2+1
for co in goodco:
    for cluster in goodclus:
        for y in prot2clus[seckeys[co]]:
            if y==cluster:
                fco.append(co)
                f3=f3+1
print(len(forth)+len(fpara)+len(fco))
#for key in indhalfkeys:
#    print(key)
#    print(sstart[key])
#    print()


#-------------------------domain code---------------------------------------

homology_search(domseq,"firstDB","domtxtHelper.txt")


helperfile=open("domtxtHelper.txt","r")
lines = helperfile.readlines()

send={}
sstart={}
for line in lines:
    b1=0
    b2=0
    b3=0
    b4=0
    se=""
    ss=""
    string=""
    for char in line:
        if b1==1:
            if char !='\t' and b2==0:
                string=string+char
            if b2 == 1:
                if char != '\t'and b3==0:
                    ss=ss+char
                if b3==1:
                    if char != '\t':
                        se=se+char
                    else:
                        if string not in sstart.keys():
                            sstart[string]=[]
                            send[string]=[]
                        sstart[string].append(ss)
                        send[string].append(se)
                        #print(string)
                        #print(ss)
                        #print(se)
                        break
                if char == '\t':
                    b3=1
            if char == '\t':
                b2=1
        if char == '\t':
            b1 =1


sseqids=[]
tax={}
fulltax={}



curfull="woop"

for line in lines:
    b=0
    b2=0
    fcur=""
    fkey=""
    ftax=""
    for char in line:
        if char=='\t':
            fkey=fcur
        if b == 0:
            if char != '\t':
                fcur=fcur+char
        else:
            if char=='#':
                b2=1

            if b2==1:
                if char=='\t':
                    break
                if char != '\t':
                    ftax=ftax+char
        if char=='\t':
            b=1
    if curfull!=fkey:
        curfull=fkey
        fulltax[fkey]=ftax





index=0

for line in lines:  #                for each line in the result in the blastp
    l=list(line)
    index=0
    tempstring=""
    helperbool=0
    b=0
    for chara in l:  #                    for each chararicter in the line
        if chara == '\t':
            for i in range(len(l)):#        record teh line until the tab char
                if l[index+i+1] != '\t':
                    tempstring= tempstring + l[index+i+1]
                else:
                    for sseqid in sseqids: #                   for each sseqid
                        if str(sseqid) == tempstring:
                            helperbool=1
                    if helperbool ==0:
                        sseqids.append(tempstring)#add ids to the array of ids
                        #print(tempstring)
                    break
            break
        index=index+1







i=0


for key in sstart.keys():
    start=(2**31-1.)
    i=i+1
    for x in sstart[key]:
        if int(x) < start:
            start=int(x)
    sstart[key]=start
    end=0
    for x in send[key]:
        if int(x) > end:
            end=int(x)
    send[key]=end
    #print(key)
    #print(sstart[key])
    #print(send[key])
print(i)

newrec=[]
for record in rec:


    for sseqid in sseqids:
        if record.id==sseqid:
            #print("-------------------------")
            #print(record.seq)
            record=SeqRecord(Seq(str(record.seq)[(sstart[sseqid]-1):(send[sseqid]+1)]), id=sseqid)

            #print("---")
            #print(record.seq)
            #print("-------------------------")
            newrec.append(record)

print(len(rec))
print(len(sseqids))
print(len(newrec))






#----------------read and writes the file for the second blastP---------------
newref=open("homologyD.fasta","w+")
SeqIO.write(newrec,"homologyD.fasta","fasta")
newref.close()

build_blast_db("homologyD.fasta","homology")

#----------------blastp of all vs all----------------------------------------
print("allvsall start")

homology_search("homologyD.fasta","homology","all_vs_allD.blastp")


print("allvsall finish")



helperfile=open("all_vs_allD.blastp","r")
lines = helperfile.readlines()

keys=[]
halfkeys=[]
indhalfkeys=[]
keysE={}
keysB={}
keysP={}
minE={}
maxB={}
maxP={}
tax={}
seckeys={}
firkeys={}
samebool=[]
taxorder={}
#sstart={}
#send={}
taxorderhelper=0

#aveP={}
#-------------go through the files and find all the info of the -------------
#--------------------AllvsAll Blastp outputfile------------------------------

#          things you should know about about blastp files
#          there differnt info is sesperated by a tab (\t)
#       in the code when you see tab(\t) its the end thing inbetween info


for line in lines:# for each line
    l=list(line)
    #print(l)
    #helper varibles
    index=0
    fin =0
    f=""
    s=""
    firststring=""
    tempstring=""
    helperbool=0
    taxbool=0
    taxa=""

    for chara in l: #                      for each character in each line
        if chara !='\t':#     if not the end of the line add to temp sting
            tempstring= tempstring +chara
            if helperbool == 0:#          find first string and add to var
                firststring= firststring+ chara
                f=f+chara
            else:
                s=s+chara
            #print(tempstring)
        if chara == '\t':#   if end of string in file add taxonomy to dic
            taxbool=0
            tax[tempstring]=taxa
            if not taxa in taxorder.keys():#       checks for no reapeats
                taxorder[taxa]=taxorderhelper
                taxorderhelper=taxorderhelper+1
            if helperbool==1:#               if the second id was passed
                if f == s: #      if first id and sencond id is the same
                    samebool.append(1)#to be honest i dont even use this
                else: #         but to scared to scared to get rid of it
                    samebool.append(0)
                index=index+1       #index helps counts the where we are
                keys.append(tempstring) #    add the full key to the dic
                charbool=0              #       We are done with id part
                secstring=""
                for char in tempstring: #   for the charartors in stirng
                    if charbool==1: #    if on second id place it in dic
                        secstring=secstring+char
                    if char=='\t':
                        charbool=1
                firkeys[tempstring]=firststring
                seckeys[tempstring]=secstring

                if not firststring in indhalfkeys:#checking for repeats
                    indhalfkeys.append(firststring)
#                      make a array of half keys to help later funtions
#                                       Both unquie and non unquie list
                halfkeys.append(firststring)
                #if not tempstring in aveB.keys():
                #    aveB[tempstring]=[]
                #    aveP[tempstring]=[]

        #------ after the ID and into the other data------
                #ss=""
                #if secstring not in sstart.keys():
                #    sstart[secstring]=[]
                for i in range(len(l)):# the start info and redoring it
                    if l[index]!='\t':
                #        ss=ss+l[index]
                        index=index+1
                    if l[index]=='\t':
                        index=index+1
                        break
                #if ss not in sstart[secstring]:
                #sstart[secstring].append(ss)
                #se=""
                #if secstring not in send.keys():
                #    send[secstring]=[]
                for i in range(len(l)):#  the end info and recording it
                    if l[index]!='\t':
                #        se=se+l[index]
                        index=index+1
                    if l[index]=='\t':
                        index=index+1
                        break
                #if se not in send[secstring]:
                #send[secstring].append(se)
                e=""#                      the E-value and recording it
                for i in range(len(l)):
                    if l[index]!='\t':
                        e=e+l[index]
                        index=index+1
                    if l[index]=='\t':
                        keysE[tempstring]=e
                        if not firststring in minE.keys():
                            minE[firststring]=e
                        else:
                            if float(e) < float(minE[firststring]):
                                minE[firststring]=e
                        index=index+1
                        break
                b=""#                   the bit-score and recording it
                for i in range(len(l)):
                    if l[index]!='\t':
                        b=b+l[index]
                        index=index+1
                    if l[index]=='\t':
                        keysB[tempstring]=b
                        if not firststring in maxB.keys():
                            maxB[firststring]=b
                        else:
                            if float(b) > float(maxB[firststring]):
                                maxB[firststring]=b

                        index=index+1
                        break
                p=""#                    The postive and recording it
                for i in range(len(l)):#   for the length of the line
                    if l[index]!='\t':#            tab inbetween info
                        p=p+l[index]
                        index=index+1 #             add more to index
                    if index==len(l):#                 if end of line
                        keysP[tempstring]=p#
                        if not firststring in maxP.keys():
                            maxP[firststring]=p#       add pos to dic
                        else:
                            if float(p) > float(maxP[firststring]):
                                maxP[firststring]=p
                        fin=1
                        # helper bool that tells us that we are done
                        index=index+1 #          gotta add the index
                        break
            tempstring=tempstring+'\t'
            helperbool=1
            taxa=taxa+'#'

        if fin == 1:#                                  if done break
            break

        if taxbool==1:     #                   stuff for the tax dic
            taxa=taxa+chara
        if chara=='#':
            taxbool =1
        index=index+1

helperfile.close()


#print(len(keys))
#for key in keys:
#    print(key)

print("allvall scanner all done")




tax2={}#                               helps dics for the taxonomy
tax1={}
taxorder={}
taxorderhelper=0

for key in keys: #                                for the all keys
    t1=""#                                               first tax
    t2=""#                                              second tax
    taxa=tax[key]#      gets the tax dic and slits it into two dic
    helperbool=0
    for char in taxa:
        if char=='#':
            helperbool=1
            continue
        if helperbool == 0:
            t1=t1+char
        if helperbool ==1:
            t2=t2+char
    if not t2 in taxorder.keys():
        taxorder[t2]=taxorderhelper
        taxorderhelper=taxorderhelper+1
    tax1[key]=t1
    tax2[key]=t2

#-----------------------finding paralogs------------------------------------
para={}
i=0
cur="woop"
helperbool=0
helperbool2=0
for key in keys:

    if cur != halfkeys[i]:
        cur=halfkeys[i]
        helperbool=0
        para[halfkeys[i]]=[]
        helperbool2=0


    if tax2[key]!=tax1[key]:
        helperbool=1
    if helperbool == 0:
        if helperbool2==1:
            #print(key)
            #print(keysE[key])
            if float(keysE[key]) <= evalue:
                para[halfkeys[i]].append(key)
            #    print(key)
    helperbool2=1
    i=i+1

#----------------------------finding orthologs---------------------------

ortho={}
i=0
cur="woop"
helperbool=0
helperbool2=0
boollist=[]
for x in range(amountofspec):
    boollist.append(0)
for key in keys:

    if cur != halfkeys[i]:
        cur=halfkeys[i]
        #helperbool=0
        ortho[halfkeys[i]]=[]
        for x in range(amountofspec):
            ortho[halfkeys[i]].append("NA")
            boollist[x]=0
        #print(key)
        #print(tax2[key])
        #print(taxorder[tax2[key]])
        boollist[taxorder[tax2[key]]]=1


    if boollist[taxorder[tax2[key]]]==0:
        if float(keysE[key]) <= evalue:
            #print(key)
            ortho[halfkeys[i]][taxorder[tax2[key]]]=key

        #    print(tax2[key])
        #    print(taxorder[tax2[key]])
            boollist[taxorder[tax2[key]]]=1

    i=i+1


#--------------------------finding co-ortholongs ---------------------------

coortho={}

for key in indhalfkeys:
    coortho[key]=[]
    for okey in ortho[key]:
        if okey!="NA":
            if seckeys[okey] in para.keys():
                for pkey in para[seckeys[okey]]:
                    coortho[key].append(pkey)
                    #print(pkey)
#-------------------spitting out the output---------------------------------


loglist=[0,0]
indexlist=[0,0]
for key in indhalfkeys:
    for orth in ortho[key]:
        if orth!="NA":
            #print(orth)
            #print(keysE[orth])
            if float(keysE[orth]) != 0.0:
                e=-(math.log10(float(keysE[orth])))

                loglist[0]=loglist[0]+e
            else:
                loglist[0]=loglist[0]+200
                e=200
            indexlist[0]=indexlist[0]+1
            #print(e)
    for corth in coortho[key]:

        #print(corth)
        if float(keysE[corth]) != 0.0:
            e=-(math.log10(float(keysE[corth])))
        #print(e)
            loglist[0]=loglist[0]+e
        else:
            loglist[0]=loglist[0]+200
            e=200
        indexlist[0]=indexlist[0]+1
    for par in para[key]:

        #print(corth)
        if float(keysE[par]) != 0.0:
            e=-(math.log10(float(keysE[par])))
        #print(e)
            loglist[1]=loglist[1]+e
        else:
            loglist[1]=loglist[1]+200
            e=200
        indexlist[1]=indexlist[1]+1
#print(loglist[0])
#print(indexlist[0])
normortho=loglist[0]/indexlist[0]
normpara=loglist[1]/indexlist[1]

output=open("graphD.abc", "w+")

network={}
for key in indhalfkeys:
    vert=[]
    for orth in ortho[key]:
        if orth!="NA":
            if float(keysE[orth]) != 0.0:
                e=-(math.log10(float(keysE[orth])))
                e=e/normortho
            else:
                e=200/normortho

            woop=str(orth+'\t'+str(e))
            output.write(woop)
            output.write("\n")
            edge=Edgy(orth,e)
            vert.append(edge)
            #print(edge.getL())
            #print(edge.getW())
    for corth in coortho[key]:
        if float(keysE[corth]) != 0.0:
            e=-(math.log10(float(keysE[corth])))
            e=e/normortho
        else:
            e=200/normortho
        woop=str(corth+'\t'+str(e))
        output.write(woop)
        output.write("\n")
        edge=Edgy(corth,e)
        vert.append(edge)
        #print(edge.getL())
        #print(edge.getW())
    for par in para[key]:
        if float(keysE[par]) != 0.0:
            e=-(math.log10(float(keysE[par])))
            e=e/normpara
        else:
            e=200/normpara
        woop=str(par+'\t'+str(e))
        output.write(woop)
        output.write("\n")
        edge=Edgy(par,e)
        vert.append(edge)
        #print(edge.getL())
        #print(edge.getW())

    network[key]=vert

output.close()

os.system("mcxload -abc graphD.abc -o OUTPUTd.mci -write-tab OUTPUTd.tab")

com = "mcl OUTPUTd.mci -I "
com = com + str(inflat)
com = com + " -o OUTPUTd.out"

os.system(com)


#---------------------
tab2prot={}
keykeys=[]

helperfile=open("OUTPUTd.tab","r")
lines = helperfile.readlines()
for line in lines:
    prot=""
    tab=""
    helperbool=0
    for char in line:
        if char == '\t':
            helperbool=1
            continue
        if helperbool==0:
            tab=tab+char
        else:
            if char=='\n':
                continue
            prot=prot+char
    #print(tab)
    #print(prot)
    tab2prot[tab]=prot
    keykeys.append(prot)


prot2cluster={}




helperfile.close()

helperfile=open("OUTPUTd.out","r")
lines = helperfile.readlines()
x=0
prot=""
cluster=""
prot2clus={}
for line in lines:
    if x<7:
        x=x+1
        continue
    prot=""
    helperbool=0
    hb=0
    if list(line)[0]!=' ':
        #prot=""
        cluster=""
        #helperbool=0
        #hb=0
        if list(line)[0] ==')':
            break
        for char in line:
            if char==" " or char=='\n':
                if helperbool==1:
                    if hb==1:
                        if prot != "$":
                            if tab2prot[prot] not in prot2clus.keys():
                                prot2clus[tab2prot[prot]]=[]
                            prot2clus[tab2prot[prot]].append(cluster)
                            #print(tab2prot[prot])
                            #print(cluster)
                prot=""
                helperbool=1
                continue
            if helperbool==0:
                cluster=cluster+char
            else:
                hb=1
                prot=prot+char
    else:
        x=0
        prot=""
        for char in line:

            if x<7:
                x=x+1
                continue

            if char==" " or char=='\n':
                if prot != "$":
                    if tab2prot[prot] not in prot2clus.keys():
                        prot2clus[tab2prot[prot]]=[]
                    prot2clus[tab2prot[prot]].append(cluster)
                    #print(tab2prot[prot])
                    #print(cluster)
                    prot=""
            else:
                prot=prot+char

goodpara=[]
goodorth=[]
goodco=[]
i1=0
i2=0
for key in indhalfkeys:

    for pkey in para[key]:
        if firkeys[pkey]!=seckeys[pkey]:
            i1=i1+1
            if firkeys[pkey] in prot2clus.keys() and seckeys[pkey] in prot2clus.keys():
                for x in prot2clus[firkeys[pkey]]:
                    for y in prot2clus[seckeys[pkey]]:
                        if x == y:
#                        print(x)
#                        print(y)
#                        print(pkey)
                            goodpara.append(pkey)
                            i2=i2+1
    for ckey in coortho[key]:
        if firkeys[ckey]!=seckeys[ckey]:
            i1=i1+1
            if firkeys[ckey] in prot2clus.keys() and seckeys[ckey] in prot2clus.keys():
                for x in prot2clus[firkeys[ckey]]:
                    for y in prot2clus[seckeys[ckey]]:
                        if x == y:
                            goodco.append(ckey)
                            i2=i2+1


    for okey in ortho[key]:
        if okey != "NA" :
            i1=i1+1
            if firkeys[okey] in prot2clus.keys() and seckeys[okey] in prot2clus.keys():
                for x in prot2clus[firkeys[okey]]:
                    for y in prot2clus[seckeys[okey]]:
                        if x == y:
                            goodorth.append(okey)
                            i2=i2+1



goodclus=[]

for key in fullrec:
    if key in prot2clus.keys():
        for clus in prot2clus[key]:
            if clus not in goodclus:

                goodclus.append(clus)

#print(goodclus)
dorth=[]
dpara=[]
dco=[]
d1=0
d2=0
d3=0
for orth in goodorth:
    for cluster in goodclus:
        for y in prot2clus[seckeys[orth]]:
            if y==cluster:
                dorth.append(orth)
                d1=d1+1
for para in goodpara:
    for cluster in goodclus:
        for y in prot2clus[seckeys[para]]:
            if y==cluster:
                dpara.append(para)
                d2=d2+1
for co in goodco:
    for cluster in goodclus:
        for y in prot2clus[seckeys[co]]:
            if y==cluster:
                dco.append(co)
                d3=d3+1

print(dorth)
print(dpara)
print(dco)
print(forth)
print(fpara)
print(fco)
print(f1)
print(f2)
print(f3)
print(d1)
print(d2)
print(d3)
'''

start=(2**31-1.)
i=0
for key in indhalfkeys:
    start=(2**31-1.)
    i=i+1
    for x in sstart[key]:
        if int(x) < start:
            start=int(x)
    sstart[key]=start
    end=0
    for x in send[key]:
        if int(x) > end:
            end=int(x)
    send[key]=end
    print(key)
    print(sstart[key])
    print(send[key])
print(i)
'''
#records = list(SeqIO.parse(goodfile, "fasta"))
#i2=0
#for x in goodpara:
    #print(x)
#    i2=i2+1
#for x in goodorth:
    #print(x)
#    i2=i2+1
#for x in goodco:
    #print(x)
#    i2=i2+1

#$print(i1)
#print(i2)
'''
if out=="":
    output=open("results.txt", "w+")
else:
    helperstring= out + ".txt"
    output=open(helperstring, "w+")
for key in indhalfkeys:
    output.write("Seq ID and File:   ")
    output.write(key)
    output.write("\n")
    output.write("PARALOGS----------------------------------------------")
    output.write("\n")
    for pkey in para[key]:
        output.write(seckeys[pkey])
        output.write("\n")
    output.write("ORTHOLOGS---------------------------------------------")
    output.write("\n")
    for okey in ortho[key]:
        if okey!="NA":
            output.write(seckeys[okey])
            output.write("\n")
    output.write("CO-ORTHOLOGS------------------------------------------")
    output.write("\n")
    for ckey in coortho[key]:
        if ckey!="NA":
            output.write(seckeys[ckey])
            output.write("\n")
    output.write("\n")
    output.write("\n")

output.close()
print(para)
#-(math.log10())

'''

'''unqcombo=[]
indexdic={}
normalize={}
for key in keys:
    if tax[key] in unqcombo:
        indexdic[tax[key]]=indexdic[tax[key]]+1
        normalize[tax[key]]=normalize[tax[key]]-(math.log10(keysE[key]))
    else:
        unqcombo.append(tax[key])
        indexdic[tax[key]]=1
        normalize[tax[key]]=(-(math.log10(keysE[key])))
'''
