import abc
from functions import *
from math import *

class CelestialObject(metaclass=abc.ABCMeta):

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        self.day_number = float(day_number)
        self.observer_latitude = float(observer_latitude)
        self.observer_longitude = float(observer_longitude)
        self.utc = float(utc)

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
    def argument_of_perihelion(self):
        pass

    @abc.abstractmethod
    def semi_major_axis(self):
        pass

    @abc.abstractmethod
    def eccentricity(self):
        pass

    @abc.abstractmethod
    def mean_anomaly(self):
        pass

    def longitude_of_perihelion(self):
        return self.longitude_of_the_ascending_node() + \
                self.argument_of_perihelion()

    def mean_longitude(self):
        return self.revolution_to_degree(self.mean_anomaly() + \
                self.longitude_of_perihelion())

    def perihelion_distance(self):
        return self.semi_major_axis() * (1 - self.eccentricity())

    def aphelion_distance(self):
        return self.semi_major_axis() * (1 + self.eccentricity())
    
    def orbital_period(self):
        return self.semi_major_axis() ** 1.5

    """ FIXME - not used unless computing comets and asteroids
    @abc.abstractmethod
    def time_of_perihelion(self):
        pass
    """

    def eccentric_anomaly(self):
        tolerance = 0.00001
        e0 = self.mean_anomaly() \
                + (180/pi) * self.eccentricity() * sind(self.mean_anomaly()) \
                * (1 + self.eccentricity() * cosd(self.mean_anomaly()))
        e1 = e0 - (e0 - (180/pi) * self.eccentricity() * sind(e0) \
                - self.mean_anomaly()) / (1 - self.eccentricity() * cosd(e0))
        while (abs(e0 - e1) > tolerance):
            e0 = e1
            e1 = e0 - (e0 - (180/pi) * self.eccentricity() * sind(e0) - \
                    self.mean_anomaly()) / (1 - self.eccentricity() * cosd(e0))
        return self.revolution_to_degree(e1)

    def x_anomaly(self):
        return self.semi_major_axis() * (cosd(self.eccentric_anomaly()) \
                - self.eccentricity())

    def y_anomaly(self):
        return self.semi_major_axis() * sind(self.eccentric_anomaly()) \
                * sqrt(1 - self.eccentricity() ** 2)

    def distance(self): 
        return sqrt(self.x_anomaly() ** 2 + self.y_anomaly() ** 2)

    def true_anomaly(self):
        return atan2d(self.y_anomaly(), self.x_anomaly())

    def revolution_to_degree(self, angle):
        return angle - floor(angle/360) * 360

    def revolution_to_hour_angle(self, angle):
        return angle - floor(angle/360 + 0.5) * 360

    def obliquity_of_the_ecliptic(self, day_number):
        return 23.4393 - 3.563e-7 * day_number

    def x_heliocentric(self):
        return self.distance() * \
                (cosd(self.longitude_of_the_ascending_node()) *\
                cosd(self.true_anomaly() + \
                self.argument_of_perihelion()) - \
                sind(self.longitude_of_the_ascending_node()) *\
                sind(self.true_anomaly() + \
                self.argument_of_perihelion()) * \
                cosd(self.inclination_to_the_ecliptic()))

    def y_heliocentric(self):
        return self.distance() * \
                (sind(self.longitude_of_the_ascending_node()) * \
                cosd(self.true_anomaly() + \
                self.argument_of_perihelion()) + \
                cosd(self.longitude_of_the_ascending_node()) * \
                sind(self.true_anomaly() + \
                self.argument_of_perihelion()) * \
                cosd(self.inclination_to_the_ecliptic()))

    def z_heliocentric(self):
        return self.distance() * \
                sind(self.true_anomaly() + \
                self.argument_of_perihelion()) * \
                sind(self.inclination_to_the_ecliptic())

    def ecliptic_longitude(self):
        return atan2d(self.y_heliocentric(), self.x_heliocentric())

    def ecliptic_latitude(self):
        return atan2d(self.z_heliocentric(), \
                sqrt(self.x_heliocentric() ** 2 + \
                self.y_heliocentric() ** 2))

    def x_heliocentric_perturbed(self):
        return self.distance() * cosd(self.ecliptic_longitude()) * \
                cosd(self.ecliptic_latitude())

    def y_heliocentric_perturbed(self):
        return self.distance() * sind(self.ecliptic_longitude()) * \
                cosd(self.ecliptic_latitude())

    def z_heliocentric_perturbed(self):
        return self.distance() * sind(self.ecliptic_latitude())

    @abc.abstractmethod
    def x_geocentric(self):
        pass

    @abc.abstractmethod
    def y_geocentric(self):
        pass

    @abc.abstractmethod
    def z_geocentric(self):
        pass

    def x_equatorial(self):
        return self.x_geocentric()

    def y_equatorial(self):
        return self.y_geocentric() * \
                cosd(self.obliquity_of_the_ecliptic(self.day_number)) -\
                self.z_geocentric() * \
                sind(self.obliquity_of_the_ecliptic(self.day_number))

    def z_equatorial(self):
        return self.y_geocentric() * \
                sind(self.obliquity_of_the_ecliptic(self.day_number)) + \
                self.z_geocentric() * \
                cosd(self.obliquity_of_the_ecliptic(self.day_number))
    
    def right_ascension(self):
        return self.revolution_to_degree(atan2d(self.y_equatorial(), \
                self.x_equatorial()))

    def declination(self):
        return atan2d(self.z_equatorial(), sqrt(self.x_equatorial() ** 2 + \
                self.y_equatorial() ** 2))

    def geocentric_distance(self):
        return sqrt(self.x_equatorial() ** 2 + self.y_equatorial() ** 2 + \
                self.z_equatorial() ** 2)

    def local_sidereal_time(self):
        return self.revolution_to_degree( \
                sun_mean_longitude(self.day_number) + 180 + \
                self.utc * 15 +  self.observer_longitude)

    # FIXME ? use revolution_to_hour_angle instead of
    # revolution_to_degree to put the hour angle between -180 and +180
    def hour_angle(self):
        return self.revolution_to_degree(self.local_sidereal_time() - \
                self.right_ascension())

    def x(self):
        return cosd(self.hour_angle()) * cosd(self.declination())

    def y(self):
        return sind(self.hour_angle()) * cosd(self.declination())

    def z(self):
        return sind(self.declination())

    def x_horizontal(self):
        return self.x() * sind(self.observer_latitude) - self.z() * \
                cosd(self.observer_latitude)

    def y_horizontal(self):
        return self.y()

    def z_horizontal(self):
        return self.x() * cosd(self.observer_latitude) + self.z() * \
                sind(self.observer_latitude)

    def azimuth(self):
        return atan2d(self.y_horizontal(), self.x_horizontal()) + 180

    def altitude(self):
        return asind(self.z_horizontal())

class Sun(CelestialObject):

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        super().__init__(day_number, observer_latitude,
                observer_longitude, utc)

    def longitude_of_the_ascending_node(self):
        return 0

    def inclination_to_the_ecliptic(self):
        return 0

    def argument_of_perihelion(self):
        return 282.9404 + 4.70935e-5 * self.day_number

    def semi_major_axis(self):
        return 1

    def eccentricity(self):
        return 0.016709 - 1.151e-9 * self.day_number

    def mean_anomaly(self): 
        return self.revolution_to_degree(356.0470 + \
                0.9856002585 * self.day_number)

    def eccentric_anomaly(self):
        return self.mean_anomaly() \
                + (180/pi) * self.eccentricity() * sind(self.mean_anomaly()) \
                * (1 + self.eccentricity() * cosd(self.mean_anomaly()))

    def x_anomaly(self):
        return cosd(self.eccentric_anomaly()) - self.eccentricity()

    def y_anomaly(self):
        return sind(self.eccentric_anomaly()) * \
                sqrt(1 - self.eccentricity() ** 2)

    """ Should work from CelestialObject class
    def distance(self): 
        return sqrt(self.x_anomaly() ** 2 + self.y_anomaly() ** 2)

    def true_anomaly(self):
        return atan2d(self.y_anomaly(), self.x_anomaly())
    """

    def longitude(self):
        return self.revolution_to_degree(self.true_anomaly() \
                + self.longitude_of_perihelion())

    def x_geocentric(self):
        return self.distance() * cosd(self.longitude())
        
    def y_geocentric(self):
        return self.distance() * sind(self.longitude())

    def z_geocentric(self):
        return 0
        
class Moon(CelestialObject):

    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        super().__init__(day_number, observer_latitude,
                observer_longitude, utc)
        self.Sun = Sun(day_number, observer_latitude, observer_longitude, utc)

    def longitude_of_the_ascending_node(self):
        return 125.1228 - 0.0529538083 * self.day_number 

    def inclination_to_the_ecliptic(self):
        return 5.1454

    def argument_of_perihelion(self):
        return 318.0634 + 0.1643573223 * self.day_number

    def semi_major_axis(self):
        return 60.2666

    def eccentricity(self):
        return 0.054900

    def mean_anomaly(self): 
        return self.revolution_to_degree(115.3654 \
                + 13.0649929509 * self.day_number)

class Mercury(CelestialObject):
    
    def __init__(self, day_number, observer_latitude,
            observer_longitude, utc):
        super().__init__(day_number, observer_latitude,
                observer_longitude, utc)
        self.Sun = Sun(day_number, observer_latitude, observer_longitude, utc)

    def longitude_of_the_ascending_node(self):
        return 48.3313 + 3.24587e-5 * self.day_number 

    def inclination_to_the_ecliptic(self):
        return 7.0047 + 5.00e-8 * self.day_number

    def argument_of_perihelion(self):
        return 29.1241 + 1.01444e-5 * self.day_number

    def semi_major_axis(self):
        return 0.387098

    def eccentricity(self):
        return 0.205635 + 5.59e-10 * self.day_number

    def mean_anomaly(self): 
        return self.revolution_to_degree(168.6562 + 4.0923344368 * \
                self.day_number)

    def x_geocentric(self):
        return self.Sun.x_geocentric() + self.x_heliocentric_perturbed() 

    def y_geocentric(self):
        return self.Sun.y_geocentric() + self.y_heliocentric_perturbed()

    def z_geocentric(self):
        return self.z_heliocentric_perturbed()
