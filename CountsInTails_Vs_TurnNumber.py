# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 13:18:20 2019

@author: msmith

Edit Log:
-------------------------------------------------------------------------------


-------------------------------------------------------------------------------



User Notes:
-------------------------------------------------------------------------------
-Code depend on a 'T###' format in the TXT file name as described below:
    -The turn number used in a text file be displated in the file name under the
        format of a capital T followed by the number wit threee numerical 
        didgets
    -Example: 'T200' and 'T003'     Bad:'t200' or 't3'
    
    -That there be NO OTHER characters using a capital T in that file name
    -Example: Good: 'Example_File_T300' and'T003 Good Name'   
               Bad: 'Example_Txt_T300' and 'This txt has T300 Turns'
    
-The txt reader skips the first 18 lines, this is becuase these are usualy text,
    not data. One can remove this function by commenting out the specified line
    within the "GetData" function

-The gaussian will only fit to the highest peak in the spectra    
    
-This code is writen assuming 1 peak per file

-Negative count numbers caused by gaussian fit overestimation are treated as 
    zeroes

-Data is not normalized, but could be if you so wished by uncommenting 4th line
    after for loop begins in Main fn
-------------------------------------------------------------------------------


"""
# imports
import numpy
import os
import math
from scipy.optimize import curve_fit
import csv
import matplotlib.pyplot as plt


# Check how many counts are in a peaks left and right tail, and prints those values to a txt file that can be read into excell for graphing



# =============================================================================
# Function definitions
# =============================================================================


#---------------------------------GetFiles-------------------------------------
# this function takes 'path', the path to the folder, and returns a list of all the txt files within

def GetFiles(path):
    files = []
    for r, d, f in os.walk(path):# r=root, d=directories, f = files
        for file in f:
            if '.txt' in file:
                files.append(os.path.join(r, file))
    return files
 
#----------------------------------GetData-------------------------------------
# this function takes a txt file path, reads all data PAST THE 18th LINE and returns the first and second colums as xData and yData respectively
          
def GetData(filePath):
    with open(filePath, newline = '') as file:
        data_reader = csv.reader(file, delimiter='\t')
        data = [line for line in data_reader]
        #line to comment out
        data = data[18:] # the first 18 lines are usualy text and not the data
        xData = [i[0] for i in data]
        yData = [i[1] for i in data]
        xData = [float(i) for i in xData]
        yData = [float(i) for i in yData]
    errorBars = numpy.asarray([math.sqrt(i+1) for i in yData])# also assigns error of yVal+1 shoudl you need that
    return [xData,yData,errorBars]

#----------------------------DataTruncation-------------------------------------
#truncate the data around the central peak for some radius along the x axis
def DataTruncation(xData,yData,radius):
    info=GetMaxInfo(xData,yData)
    yMaxIndex=info[1]
    newxData=[]
    newyData=[]
    for counter,x in enumerate(xData):
        if abs(x-xData[yMaxIndex])<radius:
            newxData.append(x)
            newyData.append(yData[counter])
    return [newxData,newyData]
              
#-----------------------------------Chi2---------------------------------------
#Calculate the chi squared value given a measurement with errors and prediction

def Chi2(y_measure,y_predict,errors):
    return numpy.sum( (y_measure - y_predict)**2 / errors**2 )

#-------------------------------Chi2reduced------------------------------------
#Calculate the reduced chi squared value given a measurement with errors and prediction, and knowing the number of parameters in the model.

def Chi2reduced(y_measure, y_predict, errors, number_of_parameters):
    return Chi2(y_measure, y_predict, errors)/(numpy.asarray(y_measure).size - number_of_parameters)

#---------------------------------LogPlot--------------------------------------
# a simple potting function that removes all values from a list below 0.5
#this is done so that a single value of 10^-10 (approx zero) gets plotted as a zero

def LogPlot(xData,yData,colour,desiredLabel,show):
    #colour and desired label are strings that specify the plot properties
    logGraphingCopy=[]
    for yValue in yData:
        if yValue<(0.5):# you can change threshold value for what constitutes 'approx zero' here
            logGraphingCopy.append(0.0)
        else:
            logGraphingCopy.append(yValue)
    
    
    plt.plot(xData,logGraphingCopy,colour,label=desiredLabel)
    plt.xlabel('Mass (u)')
    plt.ylabel('Counts')
    plt.legend()
    plt.yscale('log')
    if show:
        plt.show()

#-------------------------------Gaussfunc--------------------------------------
# a simple gaussian function
        
def Gaussfunc(x,c,sigma,a):
    return (a/(sigma*math.sqrt(numpy.pi*2)))*numpy.exp(-0.5*((x-c)/sigma)**2)

#------------------------------ReimannSum--------------------------------------
#used as integration
    
def ReimannSum(xData,yData):
    area=0
    for count in range(0,len(xData)-2):
        area+=((xData[count+1] - xData[count]) * (yData[count]))
        count+=1
    return area

#-------------------------------Normalized-------------------------------------
# normalizes the inputed data xData and yData such that the area under the curve is 1
# returns a list
    
def Normalized(xData,yData):
    normalizationFactor=ReimannSum(xData,yData)
    yData=numpy.asarray(yData)
    yData=yData/normalizationFactor
    return list(yData)

#------------------------------GetMaxInfo--------------------------------------
# Gets the index and value of the maximum yValue in a data set yData and xData
    
def GetMaxInfo(xData,yData):
    yMax=max(yData)
    yMaxIndex=yData.index(yMax)
    return [xData[yMaxIndex],yMaxIndex]

#----------------------------SubtractGaussian----------------------------------
#Simply subtracts one array from another, in this case acting as a subtration of a gaussian from our data
    
def SubtractGaussian(yData,gaussianYData):
    subtractedData=numpy.asarray(yData)-numpy.asarray(gaussianYData)
    return subtractedData

#----------------------------SubtractFWHM----------------------------------
#assumes FWHM=Sigma*2.355
#for a yData set wih a fitted gaussian of width sigma, make the y values within
#the yData set that are within a fitted gaussians FWHM range 0
    
def SubtractFWHM(xData,yData,gaussianSigma,yMaxIndex):
    FWHM=gaussianSigma*2.355
    newyData=[]
    for counter,x in enumerate(xData):
        if abs(x-xData[yMaxIndex])>(FWHM/2):
            newyData.append(yData[counter])
        else:
            newyData.append(0)
    return newyData
#----------------------------RemoveNegatives-----------------------------------
# all negative values are converted to zeroes

def RemoveNegatives(yData):
    copy=[]
    for y in yData:
        if y<0:
            copy.append(0)
        else:
            copy.append(y)
    return copy

#-------------------------------OrderFiles-------------------------------------
# this function is used if you want to prder the files by adding their index to the last character of the TXT file namr

def OrderFiles(files):
    orderedFiles=[]
    indexTracker=[]
    for currentIndex,file in enumerate(files):
        try:
            properIndex=int(file[-6:-4])-1
            indexTracker.append(properIndex)
        except:
            properIndex=int(file[-5])-1
            indexTracker.append(properIndex)
    for x in range(len(files)):
        for idx in indexTracker:
            if x==idx:
                orderedFiles.append(files[indexTracker.index(idx)])
    return orderedFiles

#-------------------------------SortList--------------------------------------
#sort a list based of ascending order of values a second reference list that is not inherently in ascending order
#assumes Max vakue in the referece list is <100000
#Both list must be the same length, and reference list CANNOT have repeat values
#both lists indecies should be coupled. ie value[4] in list1 corresponds to value[4] in list 2   
    
def SortLists(listToBeSorted, referenceList): 
    copyReferenceList=referenceList[:]# so that the reference list isnt mutated
    sortedList=[]
    sortedReferenceList=[]
    counter=len(referenceList)# so that we loop through all values in list
    while counter>0:
        #getting index of min value 
        minVal=min(copyReferenceList)
        minIndex=copyReferenceList.index(minVal)
        #assigning that min value in  to the new lists
        sortedReferenceList.append(minVal)
        sortedList.append(listToBeSorted[minIndex])
        
        #replace the min value with a high number so it does not get counted again
        copyReferenceList[minIndex]=100000
        counter-=1
    return [sortedList,sortedReferenceList]

#---------------------------GetTurnNumber--------------------------------------
#This gets the number of turns, it is dependent of T### format followed
def GetTurnNumber(file):
    file=file[-55:]
    indexT=file.find('T')
    turnNumber=int(file[indexT+1:indexT+4])
    return turnNumber

        
# =============================================================================
# Main function
# =============================================================================
    
def Main():  
    #path to the folder containing txt files
    folderPath='H:\CountsInTail\Turn Dependence test\Exported_Txt_Files'
    files=GetFiles(folderPath)#get all txt fiels in path

    leftSidePercentages=[]
    rightSidePercentages=[]
    
    #files=OrderFiles(files) # this should not be needed if T### firmat followed properly
    turnNumbers=[]
    
    for file in files:# for each individual txxt file
        print(file)
        data=GetData(file)
        xData=data[0]
        yData=data[1]
        turnNumber=GetTurnNumber(file)
        
        
        #data has to wide of a x axis, we thus remove all data outside of a 0.2 usec radius from the central peak 
        truncatedData=DataTruncation(xData,yData,0.2)
        xData=truncatedData[0]
        yData=truncatedData[1]
        #yData=Normalized(xData,yData) #Can normalize data if you want, but there is no need
        #errorBars=data[2]  #for if you want the gaussian chi squared
        integralOfData=ReimannSum(xData,yData)
        countsInData=sum(yData)
        
        yMaxInfo=GetMaxInfo(xData,yData)# list giving [LocationofMax,IndexOFMAx (in xData or yData)]
        
        #fit a gaussian to the data
        initialGuesses=[yMaxInfo[0],1,1] # syntax is [center,sigma,amplitude]
        popt, pcov = curve_fit(Gaussfunc, xData, yData, p0=initialGuesses) #returned popt is list pf best found fit perameters with same syntax as initialguesses
        gaussianFit=Gaussfunc(xData, *popt)
        
        #Find gaussian chi squared
        #NumberofPerameters=3
        #errs= numpy.sqrt(numpy.diag(pcov)) #pcov is a diagnoal matrix, to get uncertainties we take te sum of diagonal. Get a list of uncertainties with same syntax as initialguesses   
        #gaussredchi2 = Chi2reduced(yData,gaussianFit, errorBars, NumberofPerameters )
        
        #subtracting the data
        gaussianSubtractedData=SubtractGaussian(yData,gaussianFit)
        
        subtractedData=SubtractFWHM(xData,gaussianSubtractedData,popt[1],yMaxInfo[1])
        #plotting for visual understanding
        LogPlot(xData,yData,'k-','Original Data',False)
        LogPlot(xData,gaussianFit,'b-','Gaussian Fit',False)
        LogPlot(xData,subtractedData,'r-','Gaussian subtracted',True)
        
        
        #seperate the subtracted data into the left and right sides. This also removes negative counts
        leftTailxData=xData[:yMaxInfo[1]]
        leftTailyData= RemoveNegatives(subtractedData[:yMaxInfo[1]])
        rightTailxData=xData[yMaxInfo[1]:]
        rightTailyData=RemoveNegatives(subtractedData[yMaxInfo[1]:])
        
        #find number of counts in tails
        countsInLeftTail=round(sum(leftTailyData))
        countsInRightTail=round(sum(rightTailyData))
        
        #get integral of remaining counts in the left and right sides
        leftTailIntegral=ReimannSum(leftTailxData,leftTailyData)
        rightTailIntegral=ReimannSum(rightTailxData,rightTailyData)
        percentageofLeftTail=leftTailIntegral/integralOfData
        percentageofRightTail=rightTailIntegral/integralOfData
        
        
        #print statements
        print('\ntotal integral: %s'% str(integralOfData))
        print('counts in original peak: %s'% str(countsInData))
        print('----------------------')
        print('counts in left tail: %s'% str(countsInLeftTail))
        print('left tail integral: %s'% str(leftTailIntegral))
        print('left tail percentage: %s'% str(percentageofLeftTail))
        print('----------------------')
        print('counts in right tail: %s'% str(countsInRightTail))
        print('right tail integral: %s'% str(rightTailIntegral))
        print('right tail percentage: %s'% str(percentageofRightTail))
        print('--------------------------------------------------------------')
        
        
        leftSidePercentages.append(percentageofLeftTail*100)
        rightSidePercentages.append(percentageofRightTail*100)
        turnNumbers.append(turnNumber)
    
    #graphing the functions in ascending order of turn numbers
    leftSide=SortLists(leftSidePercentages,turnNumbers)
    rightSide=SortLists(rightSidePercentages,turnNumbers)
    plt.plot(leftSide[1],leftSide[0],'b--',label='Left Tail')
    plt.plot(rightSide[1], rightSide[0],'k--',label='Right Tail')
    plt.xlabel('Turn Number')
    plt.ylabel('Percentage of total area (%)')
    plt.legend()
    plt.show()    
    
    
Main()
        
        
"""
    txtFilecounter=range(len(leftSidePercentages))
    plt.plot(txtFilecounter,leftSidePercentages,'b--',label='Left Tail')
    plt.plot(txtFilecounter,rightSidePercentages,'k--',label='Right Tail')
    plt.xlabel('Turn Number')
    plt.ylabel('Pcrcentage of total area (%)')
    plt.legend()
    plt.show() 
    
    
    plt.plot(turnNumbers,leftSidePercentages,'bo',label='Left Tail')
    plt.plot(turnNumbers,rightSidePercentages,'ko',label='Right Tail')
    plt.xlabel('Turn Number')
    plt.ylabel('Pcrcentage of total area (%)')
    plt.legend()
    plt.show()
"""           
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        