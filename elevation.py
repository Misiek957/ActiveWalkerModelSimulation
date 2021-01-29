from time import sleep
import json
import urllib3
import math
import matplotlib.pyplot as plt
import cv2

select_loc = [(52.59749,-1.97889)]
coordList = []
# locations=[(52.59749,-1.97889),
#            (50.449561, 30.525366),
#            (49.449561, 31.525366)] #(lat,lon) pairs
API_KEY = 'AIzaSyBaoamIS4eiYrmP8tI9kvmtfRSE8ZXrWoQ'
loc_url = ''

## Setup Coordinates

# Assuming that the Earth is a sphere with a circumference of 40075 km.
# Length in meters of 1° of latitude = always 111.32 km
#                     0.008983° = 1km
#                     0.000 008 983 = 1m
#                     0.000 017 966 = 2m
# Length in meters of 1° of longitude = 40075 km * cos( latitude ) / 360

yLen = 111320  # length in metres latitude per degree, same for all points
yDeg = (1 / yLen)  # Degree value at that latitude indicating 2m interval
xLen = 40075000 * math.cos(select_loc[0][0]) / 360  # length in metres at latitude per degree using radian rule
xDeg = (1 / xLen)  # Degree value at longitude indicating 2m interval

# Construct an array of evenly spaced out (0.001) in a 100x100 square around the selected location
for j in range(-5, 5):
    for i in range(-5, 5):
        coX = (i*xDeg) + select_loc[0][0]
        coY = (j*yDeg) + select_loc[0][1]
        coXY = (coX, coY)
        coordList.append((coX, coY))



#%%

## Request the data

http = urllib3.PoolManager()

# for loc in locations: #Redo
loc_len = len(coordList)-1
for loc in coordList:
    loc_url = loc_url + str(loc[0]) + ',' + str(loc[1])
    if coordList.index(loc) != loc_len:
        loc_url = loc_url + '|'

try:
    # Send request
    API_URL = 'https://maps.googleapis.com/maps/api/elevation/json?locations=0,0|'+loc_url+'&key='+API_KEY
    request = http.request('GET', API_URL)  # use |
    print('Request code - ' + str(request.status))
    locData = request.data
    response = json.loads(locData)
    print('Response status - ' + str(response['status']))
    print(response)

    # Check elevation at first 3 locations
    print('Elevation at point 1 - ' + str(response['results'][0]['elevation'])) # -3492 - (0),(0)
    print('Elevation at point 2 - ' + str(response['results'][1]['elevation']))
    print('Elevation at point 3 - ' + str(response['results'][2]['elevation']))

    # Sort Elevation data into array
    areaSize = 10
    elevation = [[0 for i in range(areaSize)] for j in range(areaSize)]
    i = 0  # column
    j = 0  # row
    for k in range(len(coordList)):
        # store elevation values from the response dictionary into an array
        elevation[j][i] = response['results'][k]['elevation']
        i += 1
        if k != 0 and (k+1) % areaSize == 0:
            j += 1
            i = 0
        # plt.imshow(elevation)
        cv2.imshow('image', elevation)
        #return elevation
except ValueError:
    print("Unable to request elevation")
    # return


# print('Elevation at point 4 - ' + str(response['results'][3]['elevation']))
    # print(locData['results'])
    # response = urllib3.request.urlopen(request).read()
    # places = loads(response)
    # print('At {0} elevation is: {1}'.format(loc, places['results'][0]['elevation']))
    # sleep(1)
        # print('Error for location: {0}'.format(loc))




