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
