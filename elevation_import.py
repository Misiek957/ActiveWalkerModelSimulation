import json
import urllib3
import numpy as np
from scipy.interpolate import interp2d
import requests
import matplotlib.image as mpimg
from geopy.distance import great_circle
from geopy.distance import geodesic
import matplotlib.pyplot as plt
from scipy.interpolate import RectBivariateSpline as ReBiSpline

API_KEY = 'AIzaSyBaoamIS4eiYrmP8tI9kvmtfRSE8ZXrWoQ'
areaInterval = 12  # Creates intervals of areaInterval+1
distanceInterval = 1  # Distance in metres represented by 1 output pixel

with open("elevation_storage.json") as file:
    elevDict = json.load(file)


class Elevation:

    def __init__(self, locName, locSelectLon, locSelectLat, areaWidth, entrySelect):
        self.locName = locName
        self.locSelect = [locSelectLon, locSelectLat]  # University of Warwick - Piazza -> 52.3793,-1.5615
        self.areaWidth = areaWidth  # 420m MAX
        self.elevDictStatus = 0
        self.loc_url = ''
        self.coordList = []
        self.elevationValues = [[None] * (areaInterval+1) for i in range(areaInterval+1)]  # Setup a zero matrix
        # self.relativeGrid = np.array(np.arange())
        # Call initiation functions
        self.coord_check()
        self.entryMet = self.setup_coordinates(entrySelect)
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
            else:
                self.elevDictStatus = 0
        if self.elevDictStatus == 0:
            print("location not found in storage")  # TODO: ERROR HANDLING

    def setup_coordinates(self, entrySelect):
        # Setup Coordinates
        # TODO: Move to __init__()
        # entryMet = np.array([[0, 0]]*len(entrySelect))  # Empty array container for Entry points in M relative to center

        # yLen = 111320  # length in metres latitude per degree, same for all points
        locY = [self.locSelect[0], self.locSelect[1] + 1]
        locX = [self.locSelect[0] + 1, self.locSelect[1]]
        yLen = geodesic(self.locSelect, locY).meters  # length in metres at latitude per degree using great circle method / geodesic
        yDeg = (self.areaWidth/areaInterval)*(1 / yLen)  # Degree interval for each latitude unit square interval
        xLen = geodesic(self.locSelect, locX).meters
        # xLen = 40075000 * np.cos(np.deg2rad(self.locSelect[0])) / 360
        xDeg = (self.areaWidth/areaInterval)*(1 / xLen)  # Degree value at longitude indicating 2m interval
        self.res = [xLen, yLen]
        # Construct an array of evenly spaced out in a 100x100 square around the selected location, according to specified
        for j in range(-int(areaInterval/2), int(areaInterval/2)+1):
            for i in range(-int(areaInterval/2), int(areaInterval/2)+1):
                coX = (i*xDeg) + self.locSelect[0]
                coY = (j*yDeg) + self.locSelect[1]
                self.coordList.append((coX, coY))
        # convert Entry lon,lat points into relative values in m
        entryMet = (np.array(entrySelect) - np.array(self.locSelect))
        entryMet[:, 0] = entryMet[:, 0] * xLen
        entryMet[:, 1] = entryMet[:, 1] * yLen
        # if max value exceeds value - throw ERROR
        return entryMet

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
                    if k != 0 and ((k+1) % (areaInterval+1)==0):
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
            except IndexError:
                print('Iteration indexing error')
            else:
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

        Grid = np.arange(int(-areaInterval/2), int(areaInterval/2)+1)
        elevationInterp = interp2d(Grid, Grid, self.elevationValues, kind='cubic')  # interpolation function
        interpInterval = (distanceInterval * areaInterval) / self.areaWidth  # Scaling of the interpolation input the to extract data at exactly distanceInterval (default 1m)
        xNew = np.arange(int(-areaInterval/2), int(areaInterval/2), interpInterval)  # define the interpolation values
        yNew = np.arange(int(-areaInterval/2), int(areaInterval/2), interpInterval)
        elevationInterpNew = elevationInterp(xNew, yNew)  # store the values
        return elevationInterpNew

    def map_image_get(self):
        #     # %% Static image of current region
        locSelect = [53.061473, -4.085126]
        if 1 == 1:  # Check if file already exists
            img = mpimg.imread('saved_locations/' + str(locSelect) + '.png')
            print('location image found in storage')
        try:
            pass
        except:  # ERROR - No saved preset at location # TODO: Specify error
            API_KEY = 'AIzaSyBaoamIS4eiYrmP8tI9kvmtfRSE8ZXrWoQ'
            reqCenter = 'center=' + str(locSelect[0]) + ',' + str(locSelect[1])
            reqZoom = '&zoom=' + str(17)
            reqSize = '&size=' + str(500) + 'x' + str(500)
            reqMaptype = '&maptype=' + 'satellite'
            API_URL = 'https://maps.googleapis.com/maps/api/staticmap?' + reqCenter + reqZoom + reqSize + reqMaptype + '&key=' + API_KEY
            response = requests.get(API_URL)
            if response.ok:  # Status check
                print('Request received successfully')
            else:
                print('Request not received')
            with open('saved_locations/' + str(locSelect) + '.png', 'wb') as file:
                file.write(response.content)
            response.close()
            img = mpimg.imread('saved_locations/' + str(locSelect) + '.png')
        finally:
            plt.figure(figsize=(6, 6))
            # img = np.rot90(img)
            imgplot = plt.imshow(img)
            plt.show()
#%%

# Assuming that the Earth is a sphere with a circumference of 40075 km.
# Length in meters of 1° of latitude = always 111.32 km
#                     0.008983° = 1km
#                     0.000 008 983 = 1m
#                     0.000 017 966 = 2m
# Length in meters of 1° of longitude = 40075 km * cos( latitude ) / 360