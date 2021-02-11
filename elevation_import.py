import json
import urllib3
import math
import numpy as np
from scipy.interpolate import interp2d

API_KEY = 'AIzaSyBaoamIS4eiYrmP8tI9kvmtfRSE8ZXrWoQ'
areaInterval = 14  # only EVEN, above 14 package size error
distanceInterval = 1  # Distance in metres represented by 1 output pixel

with open("elevation_storage.json") as file:
    elevDict = json.load(file)


class Elevation:

    def __init__(self, locName, locSelectLon, locSelectLat, areaWidth):
        self.locName = locName
        self.locSelect = [locSelectLon, locSelectLat]  # University of Warwick - Piazza -> 52.3793,-1.5615
        self.areaWidth = areaWidth  # 420m MAX
        self.elevDictStatus = 0
        self.loc_url = ''
        self.coordList = []
        self.elevationValues = [[None] * areaInterval for i in range(areaInterval)]  # Setup a zero matrix
        # Call initiation functions
        self.coord_check()
        self.setup_coordinates()
        self.request_elevation()

    def coord_check(self):
        global elevDict
        # Check if coordinates have been used
        for i in range(len(elevDict['locations'])):  # TODO: check why (elevDict['locations']-1) was used before
            if (elevDict['locations'][i]['coordinate'] == self.locSelect) & \
                    (elevDict['locations'][i]['width'] == self.areaWidth) & \
                    (elevDict['locations'][i]['interval'] == areaInterval):
                print('Elevation found for ' + elevDict['locations'][i]['name'] + '\n' + str(elevDict['locations'][i]['coordinate']))
                self.elevationValues = elevDict['locations'][i]['elevation']
                self.elevDictStatus = 1
                break
        if self.elevDictStatus == 0:
            print("location not found in storage")  # TODO: ERROR HANDLING


    def setup_coordinates(self):
        # Setup Coordinates
        # TODO: Move to __init__()
        yLen = 111320  # length in metres latitude per degree, same for all points
        yDeg = (self.areaWidth/areaInterval)*(1 / yLen)  # Degree interval for each latitude unit square interval
        xLen = 40075000 * math.cos(self.locSelect[0]) / 360  # length in metres at latitude per degree using radian rule
        xDeg = (self.areaWidth/areaInterval)*(1 / xLen)  # Degree value at longitude indicating 2m interval
        self.res = [xLen, yLen]
        # Construct an array of evenly spaced out in a 100x100 square around the selected location, according to specified
        for j in range(-int(areaInterval/2), int(areaInterval/2)):
            for i in range(-int(areaInterval/2), int(areaInterval/2)):
                coX = (i*xDeg) + self.locSelect[0]
                coY = (j*yDeg) + self.locSelect[1]
                coXY = (coX, coY)
                self.coordList.append((coX, coY))

    def request_elevation(self):
        global elevDict
        if self.elevDictStatus == 0:
            ## Request the data
            http = urllib3.PoolManager()
            # elevDict = {'locations': []}  # setup new dictionary
            # for loc in locations: #Redo
            loc_len = len(self.coordList)-1
            for loc in self.coordList:
                self.loc_url = self.loc_url + str(loc[0]) + ',' + str(loc[1])
                if self.coordList.index(loc) != loc_len:
                    self.loc_url = self.loc_url + '|'
            try:
                # Send request
                API_URL = 'https://maps.googleapis.com/maps/api/elevation/json?locations=0,0|'+self.loc_url+'&key='+API_KEY
                request = http.request('GET', API_URL)  # use |
                print('Request code - ' + str(request.status))
                locData = request.data
                response = json.loads(locData)
                print('Response status - ' + str(response['status']))
                print(response)
                # Sort Elevation data into array
                # elevTemp = []
                # elevTemp = [[0 for i in range(areaInterval)] for j in range(areaInterval)]

                i = 0  # column
                j = 0  # row
                for k in range(len(self.coordList)):
                    # store elevation values from the response dictionary into an array
                    self.elevationValues[j][i] = response['results'][k+1]['elevation']
                    i += 1
                    if k != 0 and (k+1) % areaInterval == 0:
                        j += 1
                        i = 0
                # Write result to .json file storage
                elevDict['locations'].append({
                    'name': self.locName,
                    'coordinate': self.locSelect,
                    'width': self.areaWidth,
                    'resolution': self.res,
                    'interval': areaInterval,
                    'elevation': self.elevationValues
                 })
                with open('elevation_storage.json', 'w') as outfile:
                    json.dump(elevDict, outfile, indent=4)
                    outfile.write('\n')
            except ValueError:
                print("Unable to request elevation")  # TODO: ERROR HANDLING
                # return
        else:
            print("request skipped")

            # print('Elevation at point 4 - ' + str(response['results'][3]['elevation']))
                # print(locData['results'])
                # response = urllib3.request.urlopen(request).read()
                # places = loads(response)
                # print('At {0} elevation is: {1}'.format(loc, places['results'][0]['elevation']))
                # sleep(1)
                    # print('Error for location: {0}'.format(loc))

    def elev_interpolation(self):
        # Data interpolation
        xGrid = np.arange(int(-areaInterval/2), int(areaInterval/2))
        yGrid = np.arange(int(-areaInterval/2), int(areaInterval/2))
        elevationInterp = interp2d(xGrid, yGrid, self.elevationValues, kind="cubic")
        interpInterval = (distanceInterval * areaInterval) / self.areaWidth  # Scaling of the interpolation input the to extract data at exactly distanceInterval (default 1m)
        xNew = np.arange(int(-areaInterval/2), int(areaInterval/2), interpInterval)
        yNew = np.arange(int(-areaInterval/2), int(areaInterval/2), interpInterval)
        elevationInterpNew = elevationInterp(xNew, yNew)
        return elevationInterpNew
#%%

# Assuming that the Earth is a sphere with a circumference of 40075 km.
# Length in meters of 1° of latitude = always 111.32 km
#                     0.008983° = 1km
#                     0.000 008 983 = 1m
#                     0.000 017 966 = 2m
# Length in meters of 1° of longitude = 40075 km * cos( latitude ) / 360