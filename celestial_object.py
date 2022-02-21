import abc
from functions import *
from math import *

class CelestialObject:

    __metaclass__ = abc.ABCMeta

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        self.day_number = float(day_number)
        self.observer_latitude = float(observer_latitude)
        self.observer_longitude = float(observer_longitude)
        self.utc = float(utc)

    """
    def __iter__(self):
        return (i for i in (self.day_number, self.observer_latitude, 
            self.observer_longitude, self.utc))

    def __repr__(self):
        class_name = type(self).__name__
        return '{}({!r}, {!r}, {!r}, {!r})'.format(class_name, *self)
    """

    def __repr__(self):
        class_name = type(self).__name__
        return '{}({!r}, {!r}, {!r}, {!r})'.format(class_name, 
                self.day_number, self.observer_latitude,
                self.observer_longitude, self.utc)

    @abc.abstractmethod
    def longitude_of_the_ascending_node(self):
        pass
  
    @abc.abstractmethod
    def inclination_to_the_ecliptic(self):
        pass

    @abc.abstractmethod
    def argument_of_perihelion(self, day_number):
        pass

    @abc.abstractmethod
    def semi_major_axis(self):
        pass

    @abc.abstractmethod
    def eccentricity(self, day_number):
        pass

    @abc.abstractmethod
    def mean_anomaly(self, day_number):
        pass

    def longitude_of_perihelion(self):
        return self.longitude_of_the_ascending_node + \
                self.argument_of_perihelion

    def mean_longitude(self):
        return self.mean_anomaly + self.longitude_of_perihelion

    def perihelion_distance(self):
        return self.semi_major_axis * (1 - self.eccentricity)

    def aphelion_distance(self):
        return self.semi_major_axis * (1 + self.eccentricity)
    
    def orbital_period(self):
        return self.semi_major_axis ** 1.5

    @abc.abstractmethod
    def time_of_perihelion(self):
        pass

    def eccentric_anomaly(self):
        tolerance = 0.00001
        e0 = self.mean_anomaly \
                + (180/pi) * self.eccentricity * sind(self.mean_anomaly) \
                * (1 + self.eccentricity * cosd(self.mean_anomaly))
        e1 = e0 - (e0 - (180/pi) * self.eccentricity * sind(e0) \
                - self.mean_anomaly) / (1 - self.eccentricity * cosd(e0))
        while (abs(e0 - e1) > tolerance):
            e0 = e1
            e1 = e0 - (e0 - (180/pi) * self.eccentricity * sind(e0) - \
                    self.mean_anomaly) / (1 - self.eccentricity * cosd(e0))
        return self.revolution_to_degree(e1)

    def x_anomaly(self):
        return cosd(self.eccentric_anomaly()) - self.eccentricity

    def y_anomaly(self):
        return sind(self.eccentric_anomaly()) * sqrt(1 - self.eccentricity ** 2)

    def distance(self): 
        return sqrt(self.x_anomaly ** 2 + self.y_anomaly ** 2)

    def true_anomaly(self):
        return atan2d(self.y_anomaly, self.x_anomaly)

    def revolution_to_degree(self, angle):
        return angle - floor(angle/360) * 360

    def revolution_to_hour_angle(self, angle):
        return angle - floor(angle/360 + 0.5) * 360

    def obliquity_of_the_ecliptic(self, day_number):
        return 23.4393 - 3.563e-7 * day_number

class Sun(CelestialObject):

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        CelestialObject.__init__(self, day_number, observer_latitude,
                observer_longitude, utc)
        self.longitude_of_the_ascending_node = 0
        self.inclination_to_the_ecliptic = 0
        self.argument_of_perihelion = 282.9404 + 4.70935e-5 * \
                self.day_number
        self.semi_major_axis = 1
        self.eccentricity = 0.016709 - 1.151e-9 * self.day_number
        self.mean_anomaly = self.revolution_to_degree (356.0470 + \
                0.9856002585 * self.day_number)
        self.obliquity_of_the_ecliptic = \
                self.obliquity_of_the_ecliptic(self.day_number)
        self.x_anomaly = self.x_anomaly()
        self.y_anomaly = self.y_anomaly()
        self.distance = self.distance()
        self.true_anomaly = self.true_anomaly()

    def longitude(self):
        return self.revolution_to_degree(self.true_anomaly \
                + self.longitude_of_perihelion())
        
class Moon(CelestialObject):

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        CelestialObject.__init__(self, day_number, observer_latitude,
                observer_longitude, utc)
        self.longitude_of_the_ascending_node = 125.1228 \
                - 0.0529538083 * self.day_number
        self.inclination_to_the_ecliptic = 5.1454
        self.argument_of_perihelion = 318.0634 \
                + 0.1643573223 * self.day_number
        self.semi_major_axis = 60.2666
        self.eccentricity = 0.054900
        self.mean_anomaly = self.revolution_to_degree (115.3654 \
                + 13.0649929509 * self.day_number)
