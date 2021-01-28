from time import sleep
import json
import urllib3

select_loc = [(0,0)]
coordList = []
# locations=[(52.59749,-1.97889),
#            (50.449561, 30.525366),
#            (49.449561, 31.525366)] #(lat,lon) pairs
API_KEY = 'AIzaSyBaoamIS4eiYrmP8tI9kvmtfRSE8ZXrWoQ'
loc_url = ''

lonInt = 0.1 # longtitude interval - x
latInt = 0.1 # latitude interval - y
# Construct an array of evenly spaced out (0.001) in a 100x100 square around the selected location
for j in range(-5, 5):
    for i in range(-5, 5):
        coX = (i*lonInt) + select_loc[0][0]
        coY = (j*lonInt) + select_loc[0][0]
        coXY = (coX, coY)
        coordList.append((coX, coY))

        # Assuming that the Earth is a sphere with a circumference of 40075 km.
        # Length in meters of 1° of latitude = always 111.32 km
        # Length in meters of 1° of longitude = 40075 km * cos( latitude ) / 360

#%%

http = urllib3.PoolManager()

# for loc in locations: #Redo
loc_len = len(coordList)-1
for loc in coordList:
    loc_url = loc_url + str(loc[0]) + ',' + str(loc[1])
    if coordList.index(loc) != loc_len:
        loc_url = loc_url + '|'

try:
    request = http.request('GET','https://maps.googleapis.com/maps/api/elevation/json?path=0,0|'+loc_url+'&samples=3&key='+API_KEY)  # use |
    print('Request code - ' + str(request.status))
    locData = request.data
    response = json.loads(locData)
    print('Response status - ' + str(response['status']))
    print(response)
    print('Elevation at point 1 - ' + str(response['results'][0]['elevation'])) # -3492 - (0),(0)
    print('Elevation at point 2 - ' + str(response['results'][1]['elevation']))
    print('Elevation at point 3 - ' + str(response['results'][2]['elevation']))

except ValueError:
    print("Unable to request elevation at {0}".format(loc))


# print('Elevation at point 4 - ' + str(response['results'][3]['elevation']))
    # print(locData['results'])
    # response = urllib3.request.urlopen(request).read()
    # places = loads(response)
    # print('At {0} elevation is: {1}'.format(loc, places['results'][0]['elevation']))
    # sleep(1)
        # print('Error for location: {0}'.format(loc))