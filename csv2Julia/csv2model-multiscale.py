#!/usr/bin/python
import sys
import csv
import re

def csv2model(reactionfile,parameterfile,ratelawfile,outputFile,paramType="inline"):
    scanIncludesFileName="scanIncludes.jl"
    ODEDict=dict()
    print(paramType)
    if paramType == "inline":
        print('Running CSV2JuliaDiffEq with parameters hard-coded into the CSV file, \
if this is not correct, re-run with the 5th argument set to \'scan\' or \'param\'')
    elif paramType== "scan":
        print('Running CSV2JuliaDiffEq with parameters left as a function call to paramFun(n), \
for all params. We will also create a paramFun.jl file that should be included and defines all parameters. \
If this is incorrect, please re-run with 5th argument set to \'inline\' or \'param\'')
    elif paramType== "param":
        print('Running CSV2JuliaDiffEq with parameters dynamically determined by a variable, \
re-run with the 5th argument set to \'scan\' or \'inline\'')
    else:
        print('The final argument was not recognised, please choose either \'scan\', \'inline\' or \'param\'')
    print('Opening {file} as rate law file'.format(file=ratelawfile))    
    ratelaws=dict()
    delayDict=dict()
    ODEIndexDict=dict()
    #let's populate a string array of rate laws
    with open(ratelawfile,'r') as f:
        csvreader=csv.reader(f)
        #skip header row
        next(csvreader)
        for line in csvreader:
            #print('adding key:{key}, val:{val} to ratelaws dict'.format(key=line[0],val=line[0]))
            ratelaws[line[0].strip()]=line[1].strip()

    #let's populate the parameter list
    print('Opening {file} as parameters file'.format(file=parameterfile))
    parametersDict=dict()
    parametersIndexDict=dict()
    parametersIndexValueList=[]
    parametersNameList=[]
    lineIndex=1
    with open(parameterfile,'r') as f:
        csvreader=csv.reader(f)
        #skpip header row
        next(csvreader)

        for line in csvreader:
            parametersDict[line[0].strip()]=str(line[1].strip())
            parametersIndexDict[line[0].strip()]=lineIndex
            parametersIndexValueList.append(str(line[1].strip()))
            parametersNameList.append(line[0].strip())
            lineIndex=lineIndex+1


    #let's iterate through the reaction file
    print('Opening {file} as reactions file'.format(file=reactionfile))
    with open(reactionfile,'r') as f:
        csvreader=csv.reader(f)
        #skip header row
        next(csvreader)
        for line in csvreader:
            #substrate, products, kinetic law, modifiers, parameters
            substrates=line[0].strip()
            products=line[1].strip()
            kineticlaw=line[2].strip()
            modifiers=line[3].strip()
            parameters=line[4].strip()

            #split up substrates, products, modifiers and parameters by space
            substratesInThisRxn=substrates.split(' ')
            substratesInThisRxn=list(filter(None,substratesInThisRxn))
            productsInThisRxn=products.split(' ')
            productsInThisRxn=list(filter(None,productsInThisRxn))
            modifiersInThisRxn=modifiers.split(' ')
            modifiersInThisRxn=list(filter(None,modifiersInThisRxn))
            parametersInThisRxn=parameters.split(' ')
            parametersInThisRxn=list(filter(None,parametersInThisRxn))
            #print(substratesInThisRxn)


            #print(ratelaws.keys())
            #lookup kinetic law

            thisLaw=ratelaws[kineticlaw]

            #now we need to go through products substrates modifiers and variable using regular expressions
            # we will substitute in the correct values from each table
            #substrates

            splitLaw=re.split('(\[[sS]\d{0,10}\])',thisLaw)
            newLaw=[]
            substrateIndex=0
            try:
                for part in splitLaw:
                    if(part and part.startswith('[') and re.search('([sS]\d{0,10})',part)):
                        substrateIndex=int(part[2:len(part)-1])-1
                        newLaw+=substratesInThisRxn[substrateIndex]
                    else:
                        newLaw+=list(part)

                thisLaw="".join(newLaw)
                #print(newLaw)
            except:
                print('error addding substrates {substrateIndex} to reaction {line}'.format(substrateIndex=substrateIndex, line=line))

            #products
            splitLaw=re.split('(\[[pP]\d{0,10}\])',thisLaw)
            newLaw=[]
            try:
                for part in splitLaw:
                    if(part and part.startswith('[') and re.search('([pP]\d{0,10})',part)):
                        productIndex=int(part[2:len(part)-1])-1
                        newLaw+=productsInThisRxn[productIndex]
                    else:
                        newLaw+=list(part)
            except:
                print('error addding products {productIndex} to reaction {line}'.format(productIndex=productIndex, line=line))
            thisLaw="".join(newLaw)

            #modifiers
            splitLaw=re.split('(\[[mM][Oo][Dd]\d{0,10}\])',thisLaw)
            newLaw=[]
            modifierIndex=0
            try:
                for part in splitLaw:
                    if(part and part.startswith('[') and re.search('[mM][Oo][Dd]\d{0,10}',part)):
                        modifierIndex=int(part[4:len(part)-1])-1
                        thisModifier=modifiersInThisRxn[modifierIndex]
                        if(thisModifier.startswith('delay(')):
                            #cut the word delay and brackets out
                            thisModifier=thisModifier[6:len(thisModifier)-1]
                            thisModDelayProperties=thisModifier.split(',')
                            thisMod=thisModDelayProperties[0]
                            thisModDelay=thisModDelayProperties[1]
                            thisDelayIndex=str(len(delayDict))
                            newLaw+='(h(p,t-tau_'+thisMod+'_'+thisDelayIndex+')[histindex_'+thisMod+'])'
                            delayDict[thisMod+'_'+thisDelayIndex]=thisModDelay
                        else:
                            newLaw+=modifiersInThisRxn[modifierIndex]
                    else:
                        newLaw+=list(part)
            except:
                print('error addding modifiers {modifierIndex} to reaction {line}'.format(modifierIndex=modifierIndex, line=line))
            thisLaw="".join(newLaw)


            #parameters
            splitLaw=re.split('\{(\w{1,20})\}',thisLaw)
            listLength=len(splitLaw)
            newLaw=[]
            currentParamInfo=0
            try:
                for i in range(listLength):
                    if splitLaw[i]:
                        parameterAdded=0
                        for j in range(len(parametersInThisRxn)):
                            currentParamInfo=parametersInThisRxn[j]
                            thisParameterType=str.split(parametersInThisRxn[j],'_')[0]

                            #if splitLaw[i].startswith(thisParameterType):
                            if splitLaw[i]==thisParameterType:
                                if paramType=="scan":
                                    if "(t)" not in parametersDict[parametersInThisRxn[j]]:
                                        newLaw+=list("paramFun(\""+str(parametersInThisRxn[j])+"\",modify)")
                                    else:
                                        newLaw+=list(str(parametersDict[parametersInThisRxn[j]]))                                    
                                    parameterAdded=1
                                elif paramType=="param":
                                    thisParamVal=str(parametersDict[parametersInThisRxn[j]])
                                    if "(t)" in thisParamVal:
                                        newLaw+=list("p["+str(parametersIndexDict[parametersInThisRxn[j]])+"](t)")
                                    else:
                                        newLaw+=list("p["+str(parametersIndexDict[parametersInThisRxn[j]])+"]")
                                    parameterAdded=1
                                else:
                                    newLaw+=list(str(parametersDict[parametersInThisRxn[j]]))
                                    parameterAdded=1
                        else:
                            if not parameterAdded:
                                newLaw+=splitLaw[i]


            except:
                print('error addding parameters {parametersInThisRxn} to reaction {line}\n'.format(parametersInThisRxn=parametersInThisRxn, line=line))
                print('error addding parameter: {currentParamInfo}\n'.format(currentParamInfo=currentParamInfo) )
            thisLaw="".join(newLaw)


            #we need to add this reaction to every product and substrate involved in this reaction

            for thisSubstrate in substratesInThisRxn:
                if thisSubstrate in ODEDict:
                    ODEDict[thisSubstrate]=ODEDict[thisSubstrate]+' - '+thisLaw
                else:
                    ODEDict[thisSubstrate]='dy['+str(len(ODEDict)+1)+']= -'+thisLaw
                    ODEIndexDict[len(ODEDict)]=thisSubstrate
            for thisProduct in productsInThisRxn:
                if thisProduct in ODEDict:
                    ODEDict[thisProduct]=ODEDict[thisProduct]+' + '+thisLaw
                else:
                    ODEDict[thisProduct]='dy['+str(len(ODEDict)+1)+']= + '+thisLaw
                    ODEIndexDict[len(ODEDict)]=thisProduct
            #sometimes a modifier needs an ODE but has no changes other than events
            for thisModifier in modifiersInThisRxn:
                if thisModifier not in ODEDict and not thisModifier.startswith("delay("):
                #if thisModifier not in ODEDict:
                    ODEDict[thisModifier]='dy['+str(len(ODEDict)+1)+']=0'
                    ODEIndexDict[len(ODEDict)]=thisModifier
    #print(ODEDict)
    writeODEFile(ODEDict,outputFile,delayDict,ODEIndexDict,reactionfile,parameterfile,ratelawfile,len(parametersDict))
    if paramType=="scan":
        writeParamFile(scanIncludesFileName,parametersDict)
    elif paramType=="param":
        writePfile(scanIncludesFileName,parametersNameList,parametersIndexValueList)
           
def writePfile(scanIncludesFileName,parametersNameList,parametersIndexValueList):
    with open(scanIncludesFileName,'w') as f:
        f.write('#######################################################\n')
        f.write('# Generated programmatically by CSV2JuliaDiffEq.      #\n')
        f.write('# http://github.com/SiFTW/CSV2JuliaDiffEq             #\n')
        f.write('# include this file with in model running script      #\n')
        f.write('# defines all parameter values in an array p[x]       #\n')
        f.write('#      - this model style is compatable with          #\n')
        f.write('#      sensitivities etc                              #\n')
        f.write('#######################################################\n')
        f.write('\n\n')
        f.write('paramVals=[\n')
        for index,val in enumerate(parametersIndexValueList):
            f.write(str(parametersIndexValueList[index])+' #p['+str(index)+"] is "+parametersNameList[index]+"\n")        
        f.write(']')
        f.write('\n\n')
        f.write('parameterNameList=[')
        for index,val in enumerate(parametersIndexValueList):
            f.write('\"'+str(parametersNameList[index])+'\" #parameterNameList['+str(index)+"]="+str(parametersIndexValueList[index])+"\n")        
        f.write(']')
        f.write('\n\n')
        f.write('println(\"parameters can now be searched in parameterNameList by name.\")\n')
        f.write('println(\"example to modify k_binding 1.5 fold higher:\")\n')
        f.write('println(\"indexOfParam=findfirst(x->\\\"k_binding\\\"==x,parameterNameList)\")\n')
        f.write('println(\"paramVals[indexOfParam]=paramVals[indexOfParam]*1.5\")\n')        


def writeParamFile(scanIncludesFileName,parametersDict):
    with open(scanIncludesFileName,'w') as f:
        f.write('#######################################################\n')
        f.write('# Generated programmatically by CSV2JuliaDiffEq.      #\n')
        f.write('# http://github.com/SiFTW/CSV2JuliaDiffEq             #\n')
        f.write('# include this file with in model running script      #\n')
        f.write('# defines a function and a dictionary:                #\n')
        f.write('#      - \"paramFun\" which returns parameter value     #\n')
        f.write('#      - \"modify\" which can be used to only modify    #\n')
        f.write('#        certain parameters. Default is all ones.     #\n')
        f.write('# to modify param k1 by 1.5x simply do:               #\n')
        f.write('#      modify(\"paramName\")=1.5                        #\n')
        f.write('#######################################################\n')
        f.write('\n\n')
        f.write('modify=Dict(')
        for (key,val) in parametersDict.items():
            if "(t)" not in val:
                f.write('\"'+str(key)+"\"=>"+str(1.0)+", ")        
        f.write(')')
        f.write('\n\n')
        f.write('parameterList=Dict(')
        for key,val in parametersDict.items():
            if "(t)" not in val:
                f.write('\"'+str(key)+"\"=>"+str(val)+", ")        
        f.write(')')
        f.write('\n\n')
        f.write('function paramFun(paramName,modify)\n')
        f.write('   return parameterList[paramName]*modify[paramName]\n')
        f.write('end\n')
        f.write('println(\"parameters can now be modified by name.\")\n')
        f.write('println(\"example to modify k_binding 1.5 fold higher:\")\n')
        f.write('println(\"modify[\\\"k_binding\\\"]=1.5\")\n')
        


def writeODEFile(ODEDict,outputFile,delayDict,ODEIndexDict,reactionfile,parameterfile,ratelawfile,numberOfParameters,):
    #this function will write the ODE file ready to be called by Julia

    with open(outputFile,'w') as f:
        f.write('#######################################################\n')
        f.write('# Generated programmatically by CSV2JuliaDiffEq.      #\n')
        f.write('# http://github.com/SiFTW/CSV2JuliaDiffEq             #\n')
        f.write('#######################################################\n')
        f.write('# generated from:\n')
        f.write('#    reactions file: {file}\n'.format(file=reactionfile))
        f.write('#    parameters file file: {file}\n'.format(file=parameterfile))
        f.write('#    rate law file: {file}\n'.format(file=ratelawfile))
        f.write('#\n')
        f.write('# Statistics:\n')
        f.write('#    Equations:{number}\n'.format(number=len(ODEIndexDict)))
        f.write('#    Parameters:{number}\n'.format(number=numberOfParameters))
        f.write('#######################################################\n\n')
        f.write('\n\n')
        odeFileName=outputFile.split(".")[0]
        odeNameDict=dict()
        if len(delayDict)>0:
            f.write('function {name}(dy,y,h,p,t)\n'.format(name=odeFileName))
        else:
            f.write('function {name}(dy,y,p,t)\n'.format(name=odeFileName))
        #let's deal with time-dependent params
        for line in ODEIndexDict.keys():
            f.write('\t'+ODEIndexDict[line]+'=maximum([y['+str(line)+'],0])\n')
            odeNameDict[ODEIndexDict[line]]=line
        delayOdeNameList=[]
        for delayEntry in delayDict.keys():
            f.write('\ttau_'+delayEntry+'='+delayDict[delayEntry]+'\n')
            odeName=delayEntry.split('_')[0]
            delayOdeNameList.append(odeName)
        for name in delayOdeNameList:
            f.write('\thistindex_'+name+'='+str(odeNameDict[name])+'\n')
        for key in ODEDict.keys():
            f.write('\t#'+key+'\n')
            f.write('\t'+ODEDict[key]+'\n')
        f.write('end\n')

    with open('variableNames.jl','w') as f:
        f.write('syms=[')
        for line in ODEIndexDict.keys():
            f.write('\"'+ODEIndexDict[line]+'\",')
        f.write(']')

if len(sys.argv)==5:
    csv2model(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])
elif len(sys.argv)==6:
    csv2model(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])
