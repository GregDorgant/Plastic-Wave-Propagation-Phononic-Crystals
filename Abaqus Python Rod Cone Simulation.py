import time
import os
import odbAccess
from abaqusConstants import *
from numpy import (asarray, sqrt, Inf, asfarray, mean, trapz,
                   eye, zeros, squeeze)
import numpy as np
import warnings
import datetime
import mesh
import regionToolset
import visualization
from optimization import *
from job import *
import sys, math

session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)

def MaterialAndSectionDefinition(mod, rod):
    
    ## Material, Section, and Mesh
    # Create material
    Matl = rod['Material Parameters']
    A = rod['Cross-Sectional Area']
    num_uc = rod['Number of Unit Cells']
    mat = mod.Material(name = Matl['Material Name'])
    mat.Elastic(table = ((Matl['Elastic Modulus'], Matl['Poissons Ratio']), ))
    mat.Density(table = ((Matl['Density'], ), ))
    StressStrainfilename = Matl['Stress-Strain Data Filename']
    f = open(StressStrainfilename)
    D = f.read()
    f.close()
    D1 = D.split('\n')
    PlasticStressStrain = np.zeros((len(D1)-1,2))
    for i in range(0,len(D1)-1):
        temp = D1[i].split('\t')
    
        for t in range(0,2):
            PlasticStressStrain[i][t] = float(temp[t])
    PlasticStressStrain[0][1] = 0.0
    PlasticStressStrain = tuple(map(tuple,PlasticStressStrain))
    mat.Plastic(PlasticStressStrain)
    p=0
    n_el_uc = len(A)
    Rod = mod.parts['Rod']
    for j in range(0,num_uc):    
        for i in range(0,n_el_uc):
            Lab = Rod.elements[p].label
            NameIt = 'Area ' + str(i)
            El_NameIt = 'Element ' + str(Lab)
            Setty = Rod.SetFromElementLabels(name = El_NameIt, elementLabels = (Lab,))
            mod.TrussSection(name = NameIt, material = Matl['Material Name'], 
                             area = A[i])
            Rod.SectionAssignment(region = Setty,sectionName = NameIt)
            p = p+1
    Rod.SetFromNodeLabels(name = 'Measurement Nodes',nodeLabels = tuple(rod['Measurement Elements']))        
    Rod.SetFromElementLabels(name = 'Measurement Elements',elementLabels = tuple(rod['Measurement Elements']))
    return

def ApplyLoad_and_AnalysisStep(mod,rod):

    # Create step during which our load is applied

#    mod.ExplicitDynamicsStep(name = rod['Analysis Step Name'], previous = 'Initial',
#                             timePeriod = rod['Simulation Duration'],
#                             timeIncrementationMethod=FIXED_USER_DEFINED_INC,
#                             userDefinedInc= rod['Time Resolution'],
#                             linearBulkViscosity = 0.06,quadBulkViscosity = 0.12)

    mod.ExplicitDynamicsStep(name = rod['Analysis Step Name'], previous = 'Initial', timePeriod = rod['Simulation Duration'],
                             timeIncrementationMethod=AUTOMATIC_GLOBAL)

    ## Create Load    
    # Create input load data
    
    #Use the input file name specified if the amplitude profile name is not of a specific type
    
    print(rod['Input File Name'])
    f = open(rod['Input File Name'])
    D = f.read()
    f.close()
    D1 = D.split('\n')
    Force = np.zeros((len(D1)-1,2))
    for i in range(0,len(D1)-1):
        temp = D1[i].split('\t')
        
        for t in range(0,2):
            Force[i][t] = float(temp[t])
    temp = np.transpose(Force)
    Force = map(tuple,Force)
    Force = tuple(Force)
    mod.TabularAmplitude(name = 'Amp', timeSpan = STEP, data = Force)

    
    # Apply the amplitude profile specified in rod['Amplitude Profile Name'] to
    # the left end of the rod with magnitude specified in rod['Force Magnitude']
    InputPoint = regionToolset.Region(vertices = mod.rootAssembly.instances['Rod-1'].vertices[0:1])
    mod.ConcentratedForce(name = 'Load-1',region = InputPoint, 
                          createStepName = rod['Analysis Step Name'],
                          cf1 = rod['Force Magnitude']*1.00, 
                          amplitude = 'Amp')
    
    mod.FieldOutputRequest(createStepName=rod['Analysis Step Name'],
                         name='Element Measure',
                         timeInterval = EVERY_TIME_INCREMENT, variables=('S','PE','LE'),
                         region = mod.rootAssembly.allInstances['Rod-1'].sets['Measurement Elements'],
                         position=INTEGRATION_POINTS)
    return 

def SubmitJob(rod):
    job_name = rod['Job Name']

    #If you stop/abort Abaqus while it is running a job, it will keep that job
    #in a lock file. These lines will detect if that file exists then delete it
    if os.path.exists(job_name + ".lck"):
        os.remove(job_name + ".lck")
        os.remove(job_name + ".dat")
        print('Lock file detected and deleted')
    bar_job = mdb.Job(name = job_name, model = rod['Model Name'], 
                      type = ANALYSIS, numDomains = 4, numCpus = 4)
    bar_job.submit()
    bar_job.waitForCompletion()
    
    return

def barModel(rod):  
    
    mod_name = rod['Model Name']
    
    # Create model
    mod = mdb.Model(name = mod_name)       
    assem=mod.rootAssembly
    assem.DatumCsysByDefault(CARTESIAN)
        
    num_uc = rod['Number of Unit Cells'] #Number of Unit Cells in the rod
    L_uc = rod['Unit Cell Length']
    L_rod = num_uc*L_uc
    Rod = mod.Part(name = 'Rod',dimensionality = TWO_D_PLANAR,type = DEFORMABLE_BODY)
    Sketch = mod.ConstrainedSketch(name= 'Sketch-1', sheetSize=1.0)
    Sketch.Line(point1 = (0.0,0.0), point2 = (L_rod,0.0))
    Rod.BaseWire(sketch = Sketch)
    RodInst = assem.Instance(dependent=ON,name = 'Rod-1', part = Rod)

    # Seed and mesh rod
    assem.regenerate()
    A = rod['Cross-Sectional Area']
    n_el_uc = len(A)
    Rod.seedEdgeByNumber(edges = (Rod.edges[0],),number = n_el_uc*num_uc)
    elemType = mesh.ElemType(elemCode = T2D2, elemLibrary = EXPLICIT)
    Rod.setElementType(regions = (Rod.edges[0:], ), elemTypes = (elemType,))
    Rod.generateMesh()
 
    return mod

def extractOutput(rod):
    JobName = rod['Job Name']
    MeasurementPoints = rod['Measurement Elements']
    OutputFile = rod['Output File Name']
    ResStrainOutputFile = 'Residual Strain ' + rod['Output File Name']

    myOdb = openOdb(path = JobName + '.odb', readOnly=True)  
    myViewport = session.Viewport(name= 'PnC_viewport', origin=(0, 0), width=50, height=50)
    myViewport.makeCurrent()
    myViewport.maximize()
    	#### Open Odb in Viewport
    myOdb = visualization.openOdb(path= JobName + '.odb', readOnly=True)                ###WHEN READ WITH CAE
    myViewport.setValues(displayedObject=myOdb)
    for i in session.xyDataObjects.keys():
        del session.xyDataObjects[i]          
    session.xyDataListFromField(variable = (('S',INTEGRATION_POINT, ((COMPONENT, 'S11' ),)),
                                            ('PE',INTEGRATION_POINT, ((COMPONENT, 'PE11' ),)),),
                                       odb=myOdb,outputPosition = INTEGRATION_POINT, 
                                       elementLabels = (('ROD-1',tuple(MeasurementPoints),),))

    OutputData = []
    tmp = session.xyDataObjects.keys()    
    for i in range(0,len(session.xyDataObjects)):
        tmp1 = tmp[i].replace('PI: ROD-1 ','')
        tmp1 = tmp1.replace(' IP: 1','')
        tmp1 = tmp1.replace(':S11','')
        tmp1 = tmp1.replace(':PE11','')
        tmp1 = tmp1.replace(':LE11','')
        tmp1 = tmp1.replace('E: ','')
        tmp1 = tmp1.replace('N: ','')
        session.xyDataObjects.changeKey(tmp[i],tmp1)
        OutputData.append(session.xyDataObjects[tmp1])
            
    Frame = myOdb.steps[rod['Analysis Step Name']].frames[-1].fieldOutputs['PE'].values
    ElLabel = []
    ResStrainData = []
    for i in range(0,len(Frame)):
        ElLabel.append(int(Frame[i].elementLabel))   
        ResStrainData.append(Frame[i].data[0]) 
    ElLabelSorted,I = np.unique(ElLabel,return_index = True)
    ResStrainData = np.array(ResStrainData)
    ResStrainData = ResStrainData[I]
    n_Total = len(ElLabelSorted)   
    Data = tuple(map(tuple,np.concatenate((np.reshape(ElLabelSorted,(n_Total,1)),np.reshape(ResStrainData,(n_Total,1))), axis = 1)))
    ResidualStrainData = session.XYData(name = 'ResStrainAtEnd',data = Data)
    session.xyReportOptions.setValues(numberFormat = SCIENTIFIC)
    session.writeXYReport(fileName = ResStrainOutputFile,xyData = ResidualStrainData,appendMode = OFF)
    session.writeXYReport(fileName = OutputFile,xyData = OutputData,appendMode = OFF)

    Frame = myOdb.steps[rod['Analysis Step Name']].frames[-1].fieldOutputs['S'].values
    ElLabel = []
    StressData = []
    for i in range(0,len(Frame)):
        ElLabel.append(int(Frame[i].elementLabel))   
        StressData.append(Frame[i].data[0]) 
    ElLabelSorted,I = np.unique(ElLabel,return_index = True)
    StressData = np.array(StressData)
    StressData = StressData[I]
    n_Total = len(ElLabelSorted)   
    Data = tuple(map(tuple,np.concatenate((np.reshape(ElLabelSorted,(n_Total,1)),np.reshape(StressData,(n_Total,1))), axis = 1)))
    StressData1 = session.XYData(name = 'Stress Data',data = Data)
    session.xyReportOptions.setValues(numberFormat = SCIENTIFIC)

    myOdb.close()
    return 

def Wrapper(rod):
    timer = [time.time()]
    #Construct and mesh rod
    mod = barModel(rod)
    timer.append(time.time())
    
    # Create and apply the material and the position dependent sectional properties
    MaterialAndSectionDefinition(mod, rod)
    timer.append(time.time())
    
    # Create analysis step, apply load to rod, and create the data request
    ApplyLoad_and_AnalysisStep(mod,rod)
    timer.append(time.time())
    
    # Create and submit job
    SubmitJob(rod)
    timer.append(time.time())
    
    # Extract data from job
    extractOutput(rod)
    timer.append(time.time())
    
    return timer


InputType = ['Step Wave','Sine Wave'] #Input profile type
fIt =[5.0] #Input frequency
AmpIt = [5.0,350.0] #Input amplitude (stress)
Constitutive_Behavior = ['Decreasingly Hardening'] #Type of constitutive behavior

#Initialize input string and file so that we can add to it
InputString = []
InputFile = []
#For every possible combination of the above inputs, form a file name to save the run under
for p in range(0,len(Constitutive_Behavior)):
    for q in range(0,len(InputType)):        
        for i in range(0,len(AmpIt)):
            if InputType[q] == 'Sine Wave':
                for j in range(0,len(fIt)):
                    InputString.append(str(int(AmpIt[i])) + ' MPa at ' + str(fIt[j]) + ' kHz with Sine Wave Input Rod Cone Explicit Mat1.txt')
                    InputFile.append('Input ' + str(fIt[j]) + ' kHz ' + InputType[q] + '.txt')
            elif 'sine' in InputType[q]:
                for j in range(0,len(fIt)):
                    InputString.append(str(int(AmpIt[i])) + ' MPa at ' + str(fIt[j]) + ' kHz with ' + InputType[q] + ' Input ' + Constitutive_Behavior[p]+ '.txt')
                    InputFile.append('Input ' + str(fIt[j]) + ' kHz ' + InputType[q] + '.txt')
            
            elif 'Modulated Harmonic' in InputType[q]:
                for j in range(0,len(fIt)):
                    InputString.append(str(int(AmpIt[i])) + ' MPa at ' + str(fIt[j]) + ' kHz with Modulated Harmonic Wave Input Rod Cone Explicit Mat1.txt')
                    InputFile.append('Input ' + str(fIt[j]) + ' kHz ' + InputType[q] + '.txt')
            else:
                InputString.append(str(int(AmpIt[i])) + ' MPa with ' + InputType[q] + ' Input Rod Cone Explicit Mat1.txt')
                InputFile.append('Input ' + InputType[q] + '.txt')
                
#for each of the combinations, form a string that corresponds to an input curve text file   
Iterate = []
for i in range(0,len(InputString)):
    if 'kHz' not in InputString[i]:
        if InputString[i][2] == 'M':
            Iterate.append((int(InputString[i][0:2])*1e6, 1, InputFile[i], InputString[i]))
        else:
            Iterate.append((int(InputString[i][0:3])*1e6, 1, InputFile[i], InputString[i]))

    else:
        temp = InputString[i].index('kHz')
        if InputString[i][2] == 'M':
            Iterate.append((int(InputString[i][0:2])*1e6, 1e3*float(InputString[i][temp-5:temp-1]), InputFile[i], InputString[i]))
        else:
            Iterate.append((int(InputString[i][0:3])*1e6, 1e3*float(InputString[i][temp-5:temp-1]), InputFile[i], InputString[i]))
            
Iterate = tuple(Iterate)  

#Iterate through each possible run
for ij in range(0,len(Iterate)):
    st = time.time() #Begin recording time for each iteration
    E = 68.9e9 #Elastic modulus
    rho = 2700.0 #Mass density
    poiss = 0.33 #Poissons ratio
    Ce = np.sqrt(E/rho) #Elastic wave speed
    num_uc = 1 #Number of repetitions of the unit cell (only 1 for rod-cone)
    Luc = 0.25*30 #Length of unit cell (total length for this case)
    rA = 0.5 #Ratio between the max cross-sectional area and min cross-sectional area
    rL = 0.7 #Ratio between cone length and total length
    d0 = 0.167 #Diameter of maximum cross-section
    
    L1 = (1-rL)*Luc #Length of uniform section
    L2 = rL*Luc #Length of cone section
    dL = np.sqrt(rA)*d0 #Diameter of minimum cross-section   
    
    A0 = (np.pi*(d0**2))/4 #Area of maxium cross-section
    AL = (np.pi*(dL**2))/4 #Area of minimum cross-section
    
    num_nodes_uc = 80*30 #Total number of elements
    n1 = int(np.floor(rR*(1-rL)*num_nodes_uc)) #Number of elements in rod section
    n2 = int(np.floor(rL*num_nodes_uc)) #Number of elements in cone section
    
    dX = Luc/num_nodes_uc; #Spatial discretization (element length)
    ARod0 = A0*np.ones((int(n1),)) #Area profile for uniform rod
    ATaper = np.linspace(A0,AL,n2) #Area profile for cone
    A_UC = np.concatenate((ARod0,ATaper))       

    L = Luc*num_uc #Total length of system
    T = 2*L/Ce; #Simulation duration, seconds, time required for elastic wave to travel twice the system length
    
    #Non-Optimization Parameter Information to pass on to the model building function
    mat_name = 'Aluminum'
    #Determine points at which we wish to measure
    n_Total = len(A_UC)*num_uc
    MeasurementPoints = list(map(int,np.linspace(0,n_Total,61))) #61 total measurement points corresponding to every 0.125 m
    MeasurementPoints[0] = 1 
    Matl = {'Elastic Modulus':E,'Density':rho,'Material Name':mat_name,'Poissons Ratio':poiss,
            'Stress-Strain Data Filename':'C:/Users/gdorgant3/OneDrive - Georgia Institute of Technology/MATLAB/EP Waves/Stress Strain Curve Mat1.txt'}
    
    inputfilename = 'C:/Users/gdorgant3/OneDrive - Georgia Institute of Technology/MATLAB/EP Waves/Input Step Wave.txt'
    inputFileLocation = 'C:/Users/gdorgant3/OneDrive - Georgia Institute of Technology/MATLAB/EP Waves/'
    
    rod = {'Material Parameters':Matl,'Number of Unit Cells':num_uc,
           'Cross-Sectional Area':A_UC,'Unit Cell Length':Luc,
           'Job Name':'PnC','Model Name':'PnC',
           'Measurement Elements':MeasurementPoints,
           'Time Resolution':dt,'Simulation Duration':T,
           'Force Magnitude':Iterate[ij][0]*A0,'Input Frequency':Iterate[ij][1],
           'Output File Name':Iterate[ij][3],'Analysis Step Name':'Explicit_Dynamic',
           'Amplitude Profile Name':Iterate[ij][2],'Input File Name':inputFileLocation + Iterate[ij][2]}
    
    timer = Wrapper(rod)
    print('Part Creation: ' + str(round(timer[1]-timer[0],3)) + '\n')
    print('Material and Section Application: ' + str(round(timer[2]-timer[1],3)) + '\n')
    print('Analysis, Load, and Data Request: ' + str(round(timer[3]-timer[2],3)) + '\n')
    print('Job: ' + str(round(timer[4]-timer[3],3)) + '\n')
    print('Data Extraction: ' + str(round(timer[5]-timer[4],3)) + '\n')
    print('Total Elapsed Time: ' + str(round(timer[-1]-st,3)) + '\n')
    print('Iteration: ' + str(ij+1) + '/' + str(len(Iterate)))
    print('Spatial Step: ' + str(round(dX*1e3,3)) + 'E-3')
    print('Time Step: ' + str(round(dt*1e7,3)) + 'E-7')
