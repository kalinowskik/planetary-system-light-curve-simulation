from random import sample
#from tokenize import PlainToken
#from traceback import print_tb
import numpy as np
from math import pi, floor
from astropy import units as u
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
from tools import getList, GaussNoise
from astropy.visualization import quantity_support
#from __future__ import print_function
from PyAstronomy.pyasl import planck
from synphot import SpectralElement

quantity_support()





#import ipympl

def GeneratePlanetData(Planet,FailedAttemptIndexes,SAMPLE_RATE=1,ObservationPeriod=365, filter='clear'):
    '''
    Generate signal of a given planet

    Parameters
    ----------
    Planet: planet class object
        planet which signal should be simulated
    FailedAttemptIndexes: float
        a fraction of failed observations, is meant to simulate the influence of bad weather. A real number from 0 to 1. 1 means no observation took place, 0.6 means 60 % of observations failed
    SAMPLE_RATE: float
        sample rate of observations in days. 2 means we get a new data point every 2 days
    ObservationPeriod: 
        the time period of the observational campaign in days
    filter: string
        A name of one of the standard photometric filters supported by synphot or "clear". A list can be found here: https://synphot.readthedocs.io/en/latest/synphot/bandpass.html#filter
        
    Returns
    -------
    time, lumlist
        a lighcurve, time is list of times corresponding to brightness measurements listed in lumlist
    '''
    Planet.prepare_obs(filter)
    time = np.arange(0,ObservationPeriod,SAMPLE_RATE)
    phi=[]
    for t in time:
        phi.append(t/Planet.period*2*pi)
    lumlist=[]
    for j in phi:
        lum=Planet.lum_total(j)
        #lum_nounit = KeepOnlyNumber(lum)
        lumlist.append(lum.value)
    
    #taking away observations lost due to a bad weather
    time=np.delete(time, FailedAttemptIndexes)
    lumlist=np.delete(lumlist, FailedAttemptIndexes)

    return time*u.d,lumlist*u.Watt

class LC:
    def __init__(self, PlanetList, Star, SAMPLE_RATE, ObservationPeriod, name):
        self.name=name
        self.PlanetList=PlanetList
        self.Star=Star
        self.SAMPLE_RATE=SAMPLE_RATE
        self.ObservationPeriod=ObservationPeriod
        self.nObservations=floor(self.ObservationPeriod/self.SAMPLE_RATE)

    def RandomObservationIndexes(self, clouds):
        '''
        Generate random indexes of days that were cloudy 

        Parameters
        ----------
        clouds: float
            a fraction of failed observations, is meant to simulate the influence of bad weather. A real number from 0 to 1. 1 means no observation took place, 0.6 means 60 % of observations failed


        Returns
        -------
        ThoseAttemptsThatFailed:
            list with random indexes of days that were cloudy 
        '''
        nFailedObservations=floor(self.nObservations*clouds)
        sample=np.arange(0, self.nObservations, 1)
        ThoseAttemptsThatFailed=np.random.choice(sample, size=nFailedObservations, replace=False, p=None) #draw n observation attempts without return
        return ThoseAttemptsThatFailed

    def GenerateNighttimeList(self, solar_day_duration, observing_night_duration):
        '''
        Compute time periods when it was dark enough to conduct observation

        Parameters
        ----------
        solar_day_duration: float
            Solar day duration. Unit: Earth days. For example to simulate observations made on Moon it should be 29.5 because the day on Moon lasts 29.5 Earth days! For simulation of observations made from Earth just put 1
        observing_night_duration: float
            Observing night duration. Unit: Earth days. For example for observations made from the Earth equator put slighly less than 0.5

        Returns
        -------
        NighttimeList:
            A list of time periods when it was dark enough to conduct observation
        '''
        daylight_time=solar_day_duration-observing_night_duration
        nDays=floor(self.ObservationPeriod/solar_day_duration)
        NighttimeList=[]
        for i in range(nDays):
            start=solar_day_duration*i
            end=start+daylight_time
            NighttimeList.append([start,end])
        return NighttimeList


    def WithinSet(self, sets, number):
        '''
        check if a given number belongs to one of the sets

        Parameters
        ----------
        sets: list
            list of sets in a form of [[set1_start,set1_stop],[set2_start,set2_stop],...[setn_start,setn_stop]], where the values are floats
        number: float
            a number to be tested
        Returns
        -------
        True or False:
            True if the number belongs to one of the sets, False otherwise
        '''
        for set in sets:
            if set[0]<number<set[1]:
                return True
        return False

    def DaylightObservationIndexes(self, solar_day_duration, observing_night_duration):
        '''
        this function compute indexes of observations that could't happen due to daylight

        Parameters
        ----------
        solar_day_duration: float
            Solar day duration. Unit: Earth days. For example to simulate observations made on Moon it should be 29.5 because the day on Moon lasts 29.5 Earth days! For simulation of observations made from Earth just put 1
        observing_night_duration: float
            Observing night duration. Unit: Earth days. For example for observations made from the Earth equator put slighly less than 0.5
        Returns
        -------
        ImpossibleObsIndexes: list
            list with indexes of observations that could't happen due to daylight
        '''
        NighttimeList=self.GenerateNighttimeList(solar_day_duration, observing_night_duration)
        ImpossibleObsIndexes=[]
        for i in range(self.nObservations):
            time=i*self.SAMPLE_RATE
            if self.WithinSet(NighttimeList, time):
                ImpossibleObsIndexes.append(i)
        return ImpossibleObsIndexes


    def Signal(self, Clouds=0, sigma=0, filter='clear', solar_day_duration=86400, observing_night_duration=86400/2): #The Clouds parameter is a fraction of observations that we are going to loose due to 'bad weather'
        '''
        this function generates the light curve

        Parameters
        ----------
        Clouds: float
            a fraction of failed observations, is meant to simulate the influence of bad weather. A real number from 0 to 1. 1 means no observation took place, 0.6 means 60 % of observations failed
        sigma: float
            a standard deviation of gaussian noise added to the data [W]
        filter: string
            A name of one of the standard photometric filters supported by synphot or "clear". A list can be found here: https://synphot.readthedocs.io/en/latest/synphot/bandpass.html#filter
        solar_day_duration: float
            Solar day duration. Unit: Earth days. For example to simulate observations made on Moon it should be 29.5 because the day on Moon lasts 29.5 Earth days! For simulation of observations made from Earth just put 1
        observing_night_duration: float
            Observing night duration. Unit: Earth days. For example for observations made from the Earth equator put slighly less than 0.5

        Returns
        -------
        None
        '''
        
        
        FailedAttemptIndexes=self.RandomObservationIndexes(Clouds)
        FailedAttemptIndexes2=self.DaylightObservationIndexes(solar_day_duration, observing_night_duration)
        badobs=list(set(FailedAttemptIndexes).union(set(FailedAttemptIndexes2)))


        self.signals={}
        self.sigma=sigma
        tot_sig=None
        for planet in self.PlanetList:
            signal=GeneratePlanetData(planet,SAMPLE_RATE=self.SAMPLE_RATE,ObservationPeriod=self.ObservationPeriod, FailedAttemptIndexes=badobs, filter=filter)
            if tot_sig==None:
                tot_sig=signal[1]
            else:
                tot_sig=tot_sig+signal[1]
            self.signals[planet.name] = signal
        
        #adding noise to the data
        if sigma != 0:
            self.noise=GaussNoise(np.zeros(len(tot_sig)), sigma)*u.Watt
            self.totalsignal=signal[0], tot_sig+self.noise
        elif sigma==0:
            self.noise=GaussNoise(np.zeros(len(tot_sig)), 0)*u.Watt
            self.totalsignal=signal[0], tot_sig


    

    def FourierTransform(self,xupperlim=0.007, yupperlim=0.15e21):
        '''
        Performs Fast Fourier Transform and saves the plotted results as a picture

        Parameters
        ----------
        xupperlim: float
            desired plotting limit in x
        yupperlim: float
            desired plotting limit in y
        Returns
        -------
        None
        '''
        from astropy.visualization import quantity_support
        quantity_support()


        x_tot=self.totalsignal[0]
        y_tot=self.totalsignal[1]

        
        SAMPLE_RATE=1/(x_tot[len(x_tot)-1]/(len(x_tot)-1))
        yf = fft(y_tot)
        xf = fftfreq(len(y_tot), 1 / SAMPLE_RATE)

        fig, axs = plt.subplots(1)

        axs.set(xlabel='Frequency [1/d]')
        axs.plot(abs(xf)*u.d**2,abs(yf)) #*u.d**2 is to elliminate some issue with plotting, the x unit is 1/d
        axs.scatter(abs(xf)*u.d**2,abs(yf),c='magenta', marker='+')
        axs.set_xlim([0, xupperlim])
        axs.set_ylim([0, yupperlim])
        fig.savefig(fname='ft_'+str(self.name)+'_'+str(self.SAMPLE_RATE)+'_'+str(self.ObservationPeriod)+'.eps', format='eps')






    def PlotComponents(self, limit=False):
        '''
        Plot signal contribution from all the planets separatebly, show and save the figure as a picture

        Parameters
        ----------
        limit: float or False
            desired plotting limit in x, default is False which means data from entire observation campaign period will be plotted
        Returns
        -------
        None
        '''
        nPlanets=len(self.signals)
        planetlist = getList(self.signals)
        fig, axs = plt.subplots(nPlanets)
        if self.sigma!=0:
            
            
            planetlist = getList(self.signals)

        if len(self.PlanetList)>1:
        
            for ax in axs.flat:
                ax.set(xlabel='time [days]', ylabel=None)
                if limit !=False:
                    ax.set_xlim([0, limit])
                ax.label_outer()
            fig.supylabel('Bolom Lum radiated + reflected [W]')

            for i in range(nPlanets):
                x=self.signals[planetlist[i]][0]
                print(type(self.signals[planetlist[i]][1]), type(np.array(self.noise)/nPlanets))
                y=self.signals[planetlist[i]][1]+self.noise/nPlanets
                axs[i].set_title(self.PlanetList[i].name)
                axs[i].scatter(x, y, s=0.3)



        else:
            x=self.signals[planetlist[0]][0]
            y=self.signals[planetlist[0]][1]+np.array(self.noise)
            axs.set(xlabel='time [days]', ylabel='Bolom Lum radiated + reflected [W]')
            axs.set_title(self.PlanetList[0].name)
            axs.scatter(x, y, s=0.3)
        
        fig.savefig(fname='signal_'+str(self.name)+'_'+str(self.SAMPLE_RATE)+'_'+str(self.ObservationPeriod)+'.eps', format='eps')








def PerformObservationOfABlackBody(temperature, object, filter='clear'):
    '''
    Computes how much irradiance from the black body gets through the filter. also plot the black body spectrum with and without a filter, show and save the figure

    Parameters
    ----------
    temperature: float
        temperature of the black body [K]
    object: string
        name of the body for plot description
    Returns:
    -------
    Irradiance: astropy quantity
        calculated irradiance value
        
    '''
    
    #based on https://stackoverflow.com/questions/22417484/plancks-formula-for-blackbody-spectrum
    temperature=temperature.value
    # Define wavelength in meters

    lam = np.arange(1.0*1e-10, 100000.*1e-10, 1e-10)

    # Get the Planck spectrum in [W/(m**2 m)] for a given temperature
    planck_spec = planck(temperature, lam=lam)

    if filter != 'clear':
        #Get transmission curve from a filter
        band = SpectralElement.from_filter(filter) 
        lamWithUnits=lam*u.m #synphot needs units
        transmission_curve=band(lamWithUnits)
        fname=band.meta['header']['descrip']

        #perform an observation
        obs=planck_spec*transmission_curve/u.m #dont need units in this function
        integral_filter = np.sum(obs*u.m) * (lam[1] - lam[0])


        name=np.random.randint(1,100000)
        fig= 'f'+str(name)
        ax = 'a'+str(name)
        globals()[fig], globals()[ax] = plt.subplots(1)
        eval(ax).set(xlabel="Wavelength [$\AA$]", title=object)
        eval(ax).set(ylabel="Spectral radiance [$W/(m^2 m)$]")
        eval(ax).plot(lam*1e10, planck_spec, 'b-', label = 'spectrum, bol')
        eval(ax).plot(lam*1e10, obs, 'r-', label= f'spectrum, {fname} filter')
        eval(ax).legend()
        eval(fig).savefig(fname=str(object)+'_'+str(temperature)+'.eps', format='eps')
        return integral_filter*u.Watt/u.m**2

    elif filter == 'clear':
        integral_bolo = np.sum(planck_spec) * (lam[1] - lam[0])
        return integral_bolo*u.Watt/u.m**2
