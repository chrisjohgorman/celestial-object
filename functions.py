from math import *

def sind(angle):
    return sin(radians(angle))

def cosd(angle):
    return cos(radians(angle))

def tand(angle):
    return tan(radians(angle))

def asind(angle):
    return degrees(asin(angle))

def acosd(angle):
    return degrees(acos(angle))

def atand(angle):
    return degrees(atan(angle))

def atan2d(y, x):
    return degrees(atan2(y, x))

def day_number(year, month, day):
    return (367 * year 
           - (7 * (year + ((month + 9) // 12)) // 4) 
           - (3 * ((year + ((month - 9) // 7)) // 100 + 1) // 4) 
           + (275 * month // 9) + day - 730515)

def sun_mean_anomaly(day_number):
    return 356.0470 + (0.9856002585 * day_number)

def sun_argument_of_perihelion(day_number):
    return 282.9404 + (4.70935e-5 * day_number)

def sun_mean_longitude(day_number):
    return sun_mean_anomaly(day_number) + \
            sun_argument_of_perihelion(day_number)

def jupiter_mean_anomaly(day_number):
    return 19.8950 + (0.0830853001 * day_number)

def saturn_mean_anomaly(day_number):
    return 316.9670 + (0.0334442282 * day_number)

def uranus_mean_anomaly(day_number):
    return 142.5905 + (0.011725806 * day_number)

