import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpl_toolkits.basemap import Basemap

def read_profile(path_to_profile):
    profile = pd.read_csv(str(path_to_profile))
    return profile

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_latlon(array, value):
    idx = np.sum(np.abs(array-value), axis=1).argmin()
    return array[idx]

def find_len_coor(elevation_profile, lat, lon):
    #profile is a pandas dataframe of cross section elevation
    latlon_list = np.array(list(zip(elevation_profile.lat.tolist(), elevation_profile.lon.tolist())))
    latlon_found = find_nearest_latlon(latlon_list, (lat, lon))
    # lat = latlon_found[0]
    # lon = latlon_found[1]
    ind_avg = elevation_profile.loc[elevation_profile['lat']==latlon_found[0]].loc[elevation_profile['lon']==latlon_found[1]].index.tolist()[0]
    # lat_ind = elevation_profile.loc[elevation_profile['lat']==find_nearest(elevation_profile['lat'],lat)].index.tolist()[0]
    # lon_ind = elevation_profile.ix[elevation_profile['lon']==find_nearest(elevation_profile['lon'], lon)].index.tolist()[0]
    # ind_avg = int((lat_ind+lon_ind)/2)
    #print(ind_avg)
    return elevation_profile['length'][ind_avg], elevation_profile['elevation'][ind_avg]

def calc_cs_trend(profile):
    lat1 = profile['lat'][0]
    lat2 = profile['lat'][len(profile)-1]
    lon1 = profile['lon'][0]
    lon2 = profile['lon'][len(profile)-1]
    cs_trend_rad = np.arctan2(np.sin(np.deg2rad(lon2-lon1))*np.cos(np.deg2rad(lat2)),
                          np.cos(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2))-np.sin(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.cos(np.deg2rad(lon2-lon1)))
    cs_trend = np.rad2deg(cs_trend_rad)
    return cs_trend

def calc_dip_trend_diff(trend, dip_d):
    new_trend = trend%360
    if abs(new_trend-dip_d)>abs((new_trend-180)-dip_d):
        return abs((new_trend-180)-dip_d)
    else:
        return abs(new_trend-dip_d)

def calc_structures(path_to_str, profile):
    #input path to structural data CSV with columns 'lat', 'lon', 'strike' and 'dip'
    # optional 'id' column to label with site/measurement name
    structures = pd.read_csv(str(path_to_str),
                         usecols=['id', 'strike', 'dip', 'lon', 'lat'])
    structures['elevation']= pd.Series()
    structures['dip_d']= pd.Series()
    structures['length']= pd.Series()
    structures['corr_dip']= pd.Series()
    structures['slope']= pd.Series()

    for n in range(len(structures)):
        structures['length'][n], structures['elevation'][n] = find_len_coor(profile, structures['lat'][n],structures['lon'][n])
        structures['dip_d'][n] = (structures['strike'][n]+90)%360

    for n in range(len(structures)):
        gradient = np.array([np.sin(np.deg2rad(structures['dip_d'][n])),
                             np.cos(np.deg2rad(structures['dip_d'][n])),
                             np.sin(np.deg2rad(structures['dip'][n]))])
        direction = np.array([np.cos(np.deg2rad(calc_cs_trend(profile))), np.sin(np.deg2rad(calc_cs_trend(profile))), 0])
        structures['corr_dip'][n] = structures['dip'][n]*gradient.dot(direction)
        structures['slope'][n] = -abs(np.tan(np.deg2rad(structures['corr_dip'][n])))

    structures.sort_values(by='length', inplace=True)
    structures.reset_index(inplace=True,drop=True)

    return structures

def cs_figure(profile, structures, min_length=None, max_length=None,
              ymin=None, ylim = None, extrap_range = 300, contacts=None,
              locations=None, id_measurements = False, save=False, save_name='cross_section.pdf'):
    if min_length is None:
        min_length = int(min(profile.length.tolist()))
    if max_length is None:
        max_length = int(max(profile.length.tolist()))
    if ymin is None:
        ymin = int(min(profile.elevation.tolist())*0.5)
    if ylim is None:
        ylim = int(max(profile.elevation.tolist())*1.3)
    min_ind = profile.loc[profile['length']==find_nearest(np.array(profile.length.tolist()),min_length)].index.tolist()[0]
    med_ind = profile.loc[profile['length']==find_nearest(np.array(profile.length.tolist()),(max_length)/2)].index.tolist()[0]
    max_ind = profile.loc[profile['length']==find_nearest(np.array(profile.length.tolist()),max_length)].index.tolist()[0]
    min_cs_lat = profile['lat'][min_ind]
    med_cs_lat = profile['lat'][med_ind]
    max_cs_lat = profile['lat'][max_ind]
    min_cs_lon = profile['lon'][min_ind]
    med_cs_lon = profile['lon'][med_ind]
    max_cs_lon = profile['lon'][max_ind]

    if locations:
        loc_coor = {}
        for location in locations:
            loc_mid = np.average([locations[location][0], locations[location][1]])
            loc_ind = profile.loc[profile['length']==find_nearest(np.array(profile.length.tolist()),loc_mid)].index.tolist()[0]
            loc_lat = profile['lat'][loc_ind]
            loc_lon = profile['lon'][loc_ind]
            loc_coor[location] = [loc_lat, loc_lon]

    plt.figure(figsize=(12,8))
    ax1 = plt.subplot2grid((2,3), (0,0), colspan=1)
    m = Basemap(projection='merc',llcrnrlat=46.2,urcrnrlat=48,llcrnrlon=-91,
                urcrnrlon=-88, resolution='i') #lat_ts=-25
    m.drawrivers(color='#99ffff')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    parallels = [min_cs_lat]#np.arange(-90,90,2.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    meridians = [min_cs_lon]#np.arange(0.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    X, Y = m([min_cs_lon, max_cs_lon, max_cs_lon, min_cs_lon], [max_cs_lat,max_cs_lat,min_cs_lat, min_cs_lat])
    m.plot(X,Y, color='k', markersize=300)

    ax2 = plt.subplot2grid((2,3), (0,1), colspan=1)
    zoom_rescale = abs(0.5*(max(profile.length)/111100))
    m = Basemap(projection='merc',llcrnrlat=med_cs_lat-zoom_rescale,
                urcrnrlat=med_cs_lat+zoom_rescale,llcrnrlon=med_cs_lon-2*zoom_rescale,
                urcrnrlon=med_cs_lon+2*zoom_rescale, resolution='i') #lat_ts=-25
    m.drawrivers(color='#99ffff')
    m.drawcoastlines()
    m.drawmapboundary(fill_color='#99ffff')
    m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    parallels = np.arange(-90,90,2.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    meridians = np.arange(0.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    X, Y = m([profile['lon'][min_ind], profile['lon'][max_ind]],
             [profile['lat'][min_ind], profile['lat'][max_ind]])
    m.plot(X, Y, linestyle='dotted', linewidth=2, zorder=50)
    if locations:
        for location in locations:
            X, Y = m(loc_coor[location][1], loc_coor[location][0])
            m.scatter(X,Y, label=location, color='r', zorder=51)
    plt.legend()

    ax3 = plt.subplot2grid((2,3), (0,2), colspan=1)
    ax3.plot(profile.length.tolist(),
         profile.elevation.tolist())
    ax3.text(0, 1, '{:.5}\n{:.5}'.format(profile.lat.tolist()[0],
                                        profile.lon.tolist()[0]),
             transform=ax3.transAxes,
             bbox=dict(facecolor='white', edgecolor='white'), zorder=50)
    ax3.text(1, 1, '{:.5}\n{:.5}'.format(profile.lat.tolist()[len(profile)-1],
                                        profile.lon.tolist()[len(profile)-1]),
             transform=ax3.transAxes,
             bbox=dict(facecolor='white', edgecolor='white'), zorder=50)
    plt.tight_layout()

    ax4 = plt.subplot2grid((2,3), (1, 0), colspan=3)
    plt.plot(profile.length.tolist(), profile.elevation.tolist())
    contact_strat = {}
    strat=0
    structures = structures.loc[structures['length']>=min_length].loc[structures['length']<=max_length]
    structures.reset_index(inplace=True, drop=True)
    for n in range(len(structures)):
        if structures['length'][n]>=min_length and structures['length'][n]<=max_length:
            x_linspace = np.linspace(structures['length'][n], max(profile.length.tolist()))
            y = structures['slope'][n]*(x_linspace-structures['length'][n])+structures['elevation'][n]
            if n>0:
                sep=structures['length'][n]-structures['length'][n-1]
            else:
                sep=0
            plt.plot(x_linspace, y, 'r')
            if id_measurements==True:
                site_name = str(structures['id'][n])
                plt.text(x_linspace[0], y[0], site_name,
                        bbox=dict(facecolor='white', edgecolor='white'), horizontalalignment='center')
            dip_d_angle = calc_dip_trend_diff(calc_cs_trend(profile),structures['dip_d'][n])
            extrap_strat = strat
            strat += sep*np.cos(np.deg2rad(dip_d_angle))*np.sin(np.deg2rad(abs(structures['dip'][n])))
            lab_height = 30 + (n * 10)%60
            plt.text(structures['length'][n], max(y)+lab_height, '{0}'.format(int(strat)), horizontalalignment='center')
            plt.vlines(structures['length'][n], structures['elevation'][n], max(y)+lab_height-4, linestyles='dotted')
            if contacts:
                for contact in contacts:
                    if n==0:
                        contact_strat[contact] = extrap_strat+(contacts[contact]-structures['length'][n])*np.sin(np.deg2rad(abs(structures['corr_dip'][n])))
                    elif contacts[contact]<=structures['length'][n] and contacts[contact]>=structures['length'][n-1]:
                        contact_strat[contact] = extrap_strat+(contacts[contact]-structures['length'][n-1])*np.sin(np.deg2rad(abs(structures['corr_dip'][n-1])))
            if n>0 and n<len(structures)-1:
                if structures['length'][n]-structures['length'][n-1]>extrap_range:
                    for extrap in range(int(structures['length'][n-1])+extrap_range, int(structures['length'][n]), extrap_range):
                        extrap_strat += extrap_range*np.sin(np.deg2rad(abs(structures['corr_dip'][n-1])))
                        text_height = 30 + (extrap * 10)%60
                        plt.text(extrap,
                                 profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height,
                                 '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                        plt.vlines(extrap,
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height-5,
                                   linestyles='dotted')
            elif n==len(structures)-1:
                if max_length-structures['length'][n]>extrap_range:
                    extrap_strat = strat
                    for extrap in range(int(structures['length'][n])+extrap_range, max_length, extrap_range):
                        extrap_strat += extrap_range*np.sin(np.deg2rad(abs(structures['corr_dip'][n])))
                        text_height = 30 + (extrap * 10)%60
                        plt.text(extrap,
                                 profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height,
                                 '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                        plt.vlines(extrap,
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height-5,
                                   linestyles='dotted')
            elif n==0:
                if structures['length'][n]>extrap_range:
                    for extrap in range(int(min(structures.length))-extrap_range, min_length, -extrap_range):
                        extrap_strat += -extrap_range*np.sin(np.deg2rad(abs(structures['corr_dip'][0])))
                        text_height = 30 + (extrap * 10)%60
                        plt.text(extrap,
                                 profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height,
                                 '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                        plt.vlines(extrap,
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                                   profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+text_height-5,
                                   linestyles='dotted')
    if contacts:
        for contact in contacts:
            xline = float(contacts[contact])
            ymax = profile.loc[profile['length']==find_nearest(profile['length'], xline)]['elevation']
            plt.vlines(xline,0,float(ymax),colors=np.random.rand(3), linewidth=3.0, linestyles='dotted', label=contact)
            y_place = np.random.randint(ymax+0.1*(ylim-ymin),ylim-0.1*(ylim-ymin))
            plt.text(xline, y_place+5, '{0}'.format(int(contact_strat[contact])), color='g',horizontalalignment='center')
            plt.vlines(xline, ymax, y_place, linestyles='dotted')
    if locations:
        for location in locations:
            y_rng = []
            x_rng = np.linspace(locations[location][0],locations[location][1])
            for x in x_rng:
                y_rng.append(float(profile.loc[profile['length']==find_nearest(profile['length'], x)]['elevation']))
            vertices = [(x_rng[0], 0)] + list(zip(x_rng, y_rng)) + [(x_rng[-1], 0)]
            poly = Polygon(vertices, facecolor = 'r', alpha = 0.4, edgecolor='none', label=location)
            plt.gca().add_patch(poly)

    ax4.text(0, 1, '{:.5}\n{:.5}'.format(profile['lat'][min_ind],
                                         profile['lon'][min_ind]),
             color='r', transform=ax4.transAxes, bbox=dict(facecolor='white', edgecolor='white'), zorder=50)
    ax4.text(1, 1, '{:.5}\n{:.5}'.format(profile['lat'][max_ind],
                                         profile['lon'][max_ind]),
             color='r', transform=ax4.transAxes, bbox=dict(facecolor='white', edgecolor='white'), zorder=50)

    plt.xlabel('meters')
    plt.ylabel('meters')

    plt.xlim(min_length,max_length)
    plt.ylim(ymin,ylim)

    plt.title('Bearing = {:.4}°'.format(calc_cs_trend(profile)))

    if save==True:
        plt.savefig(str(save_name))

    plt.legend()
    plt.show()

def make_cross_section(profile_data, str_data, **kwargs):
    '''
    Required arguments
    =====================
    profile_data: CSV file of profile data with columns
        'length' : accumulated ground length (in meters) of profile line (ranges from 0 to profile length)
        'lat' : latitude of data point
        'lon' : longitude of data point
        'elevation' : elevation (in meters) of data point
    str_data : CSV file of structural data with columns
        'strike', 'dip', 'lon', 'lat'

    Keyword arguments
    =====================
    min_length : minimum length plotted on cross section (default=0)
    max_length : minimum length plotted on cross section (default=3000)
    ymin : minimum height plotted on cross section (default=150)
    ylim : minimum height plotted on cross section (default=300)
    extrap_range : length interval to extrapolate strat height (default=300)
    contacts : dictionary of contacts to plot
        (format is {<string of contact name>:<length at which to plot contact>})
    locations : dictionary of locations to plot
        (format is {<string of location name>:<iterable length range over which to plot location>})
    save : optional save of plot (default is False)
    save_name : name of saved plot (includes format extension; default='cross_section.pdf')
    '''
    profile = read_profile(profile_data)
    structures = calc_structures(str_data,profile)
    return cs_figure(profile, structures, **kwargs)
