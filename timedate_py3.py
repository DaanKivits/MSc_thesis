#! /usr/bin/env python
from datetime import datetime, timedelta
import calendar
from pylab import floor
from matplotlib.dates import date2num, num2date

class Mydatetime(datetime):
    """ 
    Mydatetime(year, month, day[, hour[, minute[, second[, microsecond[,tzinfo]]]]])
    
    The year, month and day arguments are required. tzinfo may be None, or an
    instance of a tzinfo subclass. The remaining arguments may be ints or longs.

    This is a datetime object, with the added methods:

    toNum(): 
        return value is a floating point number 
        which gives number of days (fraction part represents hours,
        minutes, seconds) since 0001-01-01 00:00:00 UTC
    toDec():
        return value is a floating point number 
        which gives the yearfraction (2000.0875) representation
        of the given datetime object
    """
    def toNum(self):
        """
        return value is a floating point number 
        which gives number of days (fraction part represents hours,
        minutes, seconds) since 0001-01-01 00:00:00 UTC """
        return date2num(self)
    def toDec(self):
        """ 
        return value is a floating point number 
        which gives the yearfraction (2000.0875) representation
        of the given datetime object
        """
        return date2dec(self)

def Fromdatetime(date):
    dt=date.timetuple()
    return datetime(*dt[0:6])

def increase_time(dd,**kwargs):
    """ Function increases the time by specified amount"""
    return dd+timedelta(**kwargs)

def chardate(dd,cut=8):
    return dd.strftime('%Y%m%d%H%M')[0:cut]

def timegen(sd,ed,dt):
    dd=[]
    while sd <= ed:
        dd.append(Fromdatetime(sd))
        sd=sd+dt
    return dd

def itau2datetime(itau,iyear0):
    """ Function returns a datetime object from TM5s itau times"""
    date0=datetime(iyear0,1,1,0,0,0)
    if len(itau)==1:
           itau=[itau] 
    for time in itau:
            sec=time%60
            min=(time/60)%60
            hrs=(time/3600)%24
            day=(time/86400) 
            dt=timedelta(days=day,hours=hrs,minutes=min,seconds=sec)
            yield date0+dt 

def date2dec(date):
    """ Function converts datetime object to a Decimal number time  (e.g., 1991.875 such as in IDL, CCG) """
    if not isinstance(date,list): date=[date]

    newdate=[]
    for dd in date:
        Days0=date2num(datetime(dd.year,1,1))
        if calendar.isleap(dd.year):
            DaysPerYear=366.
        else:
            DaysPerYear=365.
        DayFrac=date2num(dd)
        newdate.append(dd.year+(DayFrac-Days0)/DaysPerYear)
    if len(newdate) == 1: return newdate[0] 
    return newdate

def dec2date(dectime):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python datetime object """
    dt=num2date(dec2num(dectime)).timetuple()
    return datetime(*dt[0:7])


def dec2num(dectime):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python decimal numtime """
    if not isinstance(dectime,list): dectime=[dectime]

    newdectime=[]
    for dd in dectime:
        yr=floor(dd)
        Days0=date2num(datetime(int(yr),1,1))
        if calendar.isleap(yr):
            DaysPerYear=366.
        else:
            DaysPerYear=365.
        DayFrac=(dd-yr)*DaysPerYear
        newdectime.append(Days0+DayFrac)
    if len(newdectime) == 1: return newdectime[0] 
    return newdectime

def num2dec(numtime):
    """ Function converts python decimal numtime to an IDL decimal time """
    res=num2mydate(numtime).toDec() 
    return res



def num2mydate(num):
    """ Function converts decimal time from year fraction (e.g., 1991.875 such as in IDL, CCG) to a python datetime object """
    dt=num2date(num).timetuple()
    return datetime(*dt[0:7])

def monthgen(sd,ed):
    """ Generate sequence of datetime objects spaced by one month"""
    from pylab import arange
    if ed<sd: 
        raise ValueError('start date exceeds end date')
        sys.exit(2)
    dates=[]
    for year in arange(sd.year,ed.year+2):
        for month in arange(1,13):
            date=datetime(year,month,1)
            if date > ed: return dates
            else: dates.append(date)

def nextmonth(dd):
    """ Find next 1st of the month following the date dd"""

    if dd.month == 12: 
        cc=dd.replace(year=dd.year+1)
        ee=cc.replace(month=1)
    else:
        ee=dd.replace(month=dd.month+1)
    ff = ee.replace(day=1)
    return ff


if __name__ == '__main__':
    #print monthgen(datetime(2000,1,1),datetime(2006,5,1))
    dd=datetime(2002,3,1)
    print(nextmonth(dd),dd)

        
         
