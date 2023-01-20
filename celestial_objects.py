from email import header
from math import pi, sin
from constants import *
from lightcurve import GeneratePlanetData, PerformObservationOfABlackBody
from pandas import DataFrame



class star:
    '''
    a class storing informations about a given star
    Parameters
    ----------
    luminosity: float
        star's luminosity in solar luminosities
    temperature: float
        star's photosphere temperature in Kelvins
    radius: float
        star's radius, unit: meters

    '''
    def __init__(self, luminosity, temperature, radius):
        self.luminosity = luminosity*sun_lum
        self.temperature = temperature * Kelvin    
        self.radius=radius * Meter

class planet:
    '''
    a class modelling a planet following Charpinet et al 2011, Supplementary Information, section D (page 18) (https://www.nature.com/articles/nature10631#Sec2)

    Parameters
    ----------
    albedo: float
        albedo of the planet, a number from 0 to 1
    radius: float
        radius of the planet in Earth radii
    distance: float
        a radius of the planet's orbit in AU
    beta: float
        measure of heat redistribution of the planet, where 1 means that the temperature is the same on both sides of the planet and 0 means that all the heat is kept on the illuminated side
    star: star class object
    period: float
        orbital period of the planet in Earth days
    name: string
        a desired name of the planet

    '''
    def __init__(self, albedo, radius, distance, beta, i, Star, period, name):

        self.name=name
        self.albedo = albedo
        self.distance = distance * au
        self.radius = radius * r_earth
        self.beta = beta
        self.Star=Star
        self.period = period
        self.inclination = i
        self.lightsource = Star.luminosity
        self.star_constant=self.lightsource / ( 4 * pi * (self.distance) ** 2) #Sun constant analog
        self.hot_side_temperature = ((1 - albedo) / (8 * pi * sigma) \
            * (self.lightsource / self.distance ** 2) * (1 / (1 + beta ** 4))) ** (1/4)
        self.dark_side_temperature = beta * self.hot_side_temperature
        message=DataFrame(data=[['sun constant equivalent','hot side','dark side'],[self.star_constant, self.hot_side_temperature, self.dark_side_temperature]])
        print(message)

    def area_illuminated(self, phi):
        return pi * self.radius ** 2 / 2 * (1+sin( self.inclination) * sin(phi))
    
    def area_dark(self, phi):
        return pi*self.radius**2-self.area_illuminated(phi)

    def prepare_obs(self, filter):
        '''
        This function computes the irradiance of both sides of the planet
        '''
        self.power_of_the_bright_side=PerformObservationOfABlackBody(temperature=self.hot_side_temperature, object=self.name+', bright side', filter=filter)
        self.power_of_the_dark_side=PerformObservationOfABlackBody(temperature=self.dark_side_temperature, object=self.name+', dark side', filter=filter)
        self.star_constant_withfilter=PerformObservationOfABlackBody(temperature=self.Star.temperature, object='Star', filter=filter)*(self.Star.radius/self.distance)**2

    #LUMINOSITY RADIATED
    def lum_rad(self, phi):
        return 4*(self.area_illuminated(phi)*self.power_of_the_bright_side+self.area_dark(phi)*self.power_of_the_dark_side)

    #LUMINOSITY REFLECTED
    def lum_ref(self, phi):
        return  self.area_illuminated(phi) * self.albedo * self.star_constant_withfilter

    #TOTAL LUMINOSITY
    def lum_total(self,phi):
        return self.lum_rad(phi)+self.lum_ref(phi)