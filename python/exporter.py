#----------------------------------------------------------
# File meshes.py
#----------------------------------------------------------
from progressbar import *
import os
import struct
import sys
from math import radians
 
def getVerticesFromFile( filename ):
    file_object = open(filename,'r')
    vertices = list()
    for line in file_object :
        data = line.rstrip().split(' ') 
        vertex = map(float, data)
        vertices.append(vertex)
    return vertices

def writeBphysHeader( baseFileName, particleNumber, beginFrame, endFrame, lifeTime) :
    #Build a correct filename
    paddingSize = len(str(particleNumber))
    index = '0'.zfill(max(6,paddingSize))
    stream = '00'
    fileName = baseFileName + '_' + index + '_' + stream + '.bphys'

    #Create the file binary mode
    file = open(fileName, "wb")
    file.write( struct.pack('8s', 'BPHYSICS') )
    #Compression flag
    file.write( struct.pack('I', 1) )
    #Particle number
    file.write( struct.pack('I', particleNumber) )
    #Unknown
    file.write( struct.pack('I', 64) )
    #begin, end, life
    for i in range(particleNumber) :
        file.write( struct.pack('3f', beginFrame, endFrame,  lifeTime ) )
    return 0

def exportToBphys( baseFileName, frameNumber, position, velocity) :
    
    #Build a correct filename
    particleNumber = len(position)
    paddingSize = len(str(particleNumber))
    index = str(frameNumber).zfill(max(6,paddingSize))
    stream = '00'
    fileName = baseFileName + '_' + index + '_' + stream + '.bphys'

    #Create the file binary mode
    file = open(fileName, "wb")
    file.write( struct.pack('8s', 'BPHYSICS') )
    #Compression flag
    file.write( struct.pack('I', 1) )
    #Particle number
    file.write( struct.pack('I', particleNumber) )
    #Data information : index + position + velocity = 7
    file.write( struct.pack('I', 7) )
    for i in range(particleNumber) :
        file.write( struct.pack('I', i) )
        file.write( struct.pack('f', position[i][0] ) )
        file.write( struct.pack('f', position[i][2] ) )
        file.write( struct.pack('f', position[i][1] ) )
        file.write( struct.pack('f', velocity[i][0]) )
        file.write( struct.pack('f', velocity[i][1]) )
        file.write( struct.pack('f', velocity[i][2]) )
    return 0

if __name__ == "__main__":
    #Data location
    positionDir = "../build/output/position/"
    positionData = sorted(os.listdir(positionDir))
    velocityDir = "../build/output/velocity/"
    velocityData = sorted(os.listdir(velocityDir))
    #Bphys param
    baseFileName = "../build/output/bphys/particle"
    beginFrame = 1
    lifeTime = len(velocityData)
    endFrame = beginFrame+lifeTime
    
    sample = getVerticesFromFile( positionDir+positionData[0] )
    particleNumber = len(sample)

    widgets = ['Converting : ', Percentage(), ' ', Bar(marker='+',left='[',right=']'),
               ' ', ETA(), ' ', FileTransferSpeed()] #see docs for other options
    pbar = ProgressBar(widgets=widgets, maxval=len(positionData))
    pbar.start()

    writeBphysHeader( baseFileName, particleNumber, beginFrame, endFrame, lifeTime)
    for i in range(len(positionData)) :
        position = getVerticesFromFile( positionDir+positionData[i] )
        velocity = getVerticesFromFile( velocityDir+velocityData[i] )
        particleNumber = len(position)
        exportToBphys(baseFileName, i+1, position, velocity)
        pbar.update(i) #this adds a little symbol at each iteration
    pbar.finish()
    print

