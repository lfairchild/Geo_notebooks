import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from matplotlib.patches import PathPatch
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from mpl_toolkits.basemap import Basemap

def read_profile(path_to_profile):
    profile = pd.read_csv(str(path_to_profile))
    return profile

def find_nearest(array,value):
    idx = (np.abs(array-value)).argmin()
    return array[idx]

def find_nearest_latlon(array, values):
    lat1 = values[0]
    c=[]
    for z in array:
        lat_diff = np.abs(z[0]-values[0])
        lon_diff = np.abs(z[1]-values[1])
        lat2 = z[0]
        a = (np.sin(np.deg2rad(lat_diff)/2)**2)+np.cos(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*(np.sin(np.deg2rad(lon_diff)/2)**2)
        c.append(2*np.arctan2(np.sqrt(a), np.sqrt(1-a)))
    # c is angular distance in radians
    idx = np.array(c).argmin()
    # convert c to distance in meters
    distance = np.array(c)[idx]*6371000
    # idx = np.sum(np.abs(array-value), axis=1).argmin()
    return array[idx], distance

def find_len_coor(elevation_profile, lat, lon, dip_d):
    #profile is a pandas dataframe of cross section elevation
    profile_trend = calc_cs_trend(elevation_profile)
    dip_trend_diff = calc_dip_trend_diff(profile_trend, dip_d)
    latlon_list = np.array(list(zip(elevation_profile.lat.tolist(), elevation_profile.lon.tolist())))
    nearest_latlon, smallest_distance = find_nearest_latlon(latlon_list, (lat, lon))
    projection_trend_list = np.array([(lat, lon), (nearest_latlon[0], nearest_latlon[1])])
    projection_trend = calc_trend(projection_trend_list)
    # lat = latlon_found[0]
    # lon = latlon_found[1]
    ind_avg = elevation_profile.loc[elevation_profile['lat']==nearest_latlon[0]].loc[elevation_profile['lon']==nearest_latlon[1]].index.tolist()[0]
    # lat_ind = elevation_profile.loc[elevation_profile['lat']==find_nearest(elevation_profile['lat'],lat)].index.tolist()[0]
    # lon_ind = elevation_profile.ix[elevation_profile['lon']==find_nearest(elevation_profile['lon'], lon)].index.tolist()[0]
    # ind_avg = int((lat_ind+lon_ind)/2)
    projection_trend_diff = abs(dip_d-projection_trend)
    if projection_trend_diff > 90:
        new_length = elevation_profile['length'][ind_avg]+smallest_distance*abs(np.tan(np.deg2rad(dip_trend_diff)))
        new_ind = elevation_profile.loc[elevation_profile['length']==find_nearest(elevation_profile['length'], new_length)].index.tolist()[0]
        return new_length, elevation_profile['elevation'][new_ind]
    elif projection_trend_diff < 90:
        new_length = elevation_profile['length'][ind_avg]-smallest_distance*abs(np.tan(np.deg2rad(dip_trend_diff)))
        new_ind = elevation_profile.loc[elevation_profile['length']==find_nearest(elevation_profile['length'], new_length)].index.tolist()[0]
        return new_length, elevation_profile['elevation'][new_ind]
    else:
        return elevation_profile['length'][ind_avg], elevation_profile['elevation'][ind_avg]

def calc_trend(latlon):
    lat1 = latlon[0][0]
    lat2 = latlon[1][0]
    lon1 = latlon[0][1]
    lon2 = latlon[1][1]
    trend_rad = np.arctan2(np.sin(np.deg2rad(lon2-lon1))*np.cos(np.deg2rad(lat2)),
                          np.cos(np.deg2rad(lat1))*np.sin(np.deg2rad(lat2))-np.sin(np.deg2rad(lat1))*np.cos(np.deg2rad(lat2))*np.cos(np.deg2rad(lon2-lon1)))
    trend = np.rad2deg(trend_rad)
    return trend


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
    diff1 = abs(new_trend-dip_d)
    diff2 = abs((new_trend-180)-dip_d)
    if diff1>diff2:
        return diff2
    else:
        return diff1

def calc_structures(path_to_str, profile):
    #input path to structural data CSV with columns 'lat', 'lon', 'strike' and 'dip'
    # optional 'id' column to label with site/measurement name
    try:
        structures = pd.read_csv(str(path_to_str),
                             usecols=['id', 'strike', 'dip', 'lon', 'lat'])
    except:
        structures = pd.read_csv(str(path_to_str),
                             usecols=['strike', 'dip', 'lon', 'lat'])
    structures['elevation']= pd.Series()
    structures['dip_d']= pd.Series()
    structures['length']= pd.Series()
    structures['corr_dip']= pd.Series()
    structures['slope']= pd.Series()

    for n in range(len(structures)):
        structures['dip_d'][n] = (structures['strike'][n]+90)%360
        structures['length'][n], structures['elevation'][n] = find_len_coor(profile, structures['lat'][n],structures['lon'][n], structures['dip_d'][n])

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

def calc_strat(length, trend_diff_deg, dip_degrees, elevation_diff):
    trend_diff = np.deg2rad(trend_diff_deg)
    dip = np.deg2rad(abs(dip_degrees))
    strat_height = (length*np.cos(trend_diff))*np.sin(dip)
    elevation_correction = elevation_diff*np.cos(dip)
    return strat_height+elevation_correction

def cs_figure(profile, structures, min_length=None, max_length=None,
              ymin=None, ylim = None, extrap_range = None, contacts=None,
              locations=None, id_measurements = False, reset_strat0=None, save=False, save_name='cross_section.pdf', return_all_strat_heights = False):
    all_strat_heights = {}
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

    if extrap_range is None:
        extrap_range = (max_length-min_length)//10

    profile_trend = calc_cs_trend(profile)

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
                urcrnrlon=-88, resolution='l') #lat_ts=-25
    # m.drawrivers(color='#99ffff')
    # m.drawcoastlines()
    m.drawmapboundary(fill_color='#99ffff')
    m.readshapefile('./shapefiles/MI_outline', 'MI_outline')
    m.readshapefile('./shapefiles/WI_outline', 'WI_outline')
    m.readshapefile('./shapefiles/MN_outline', 'MN_outline')
    patches   = []
    for info, shape in zip(m.MI_outline_info, m.MI_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax1.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))
    patches   = []
    for info, shape in zip(m.WI_outline_info, m.WI_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax1.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))
    patches   = []
    for info, shape in zip(m.MN_outline_info, m.MN_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax1.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))
    # m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    parallels = [min_cs_lat]#np.arange(-90,90,2.)
    m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    meridians = [min_cs_lon]#np.arange(0.,360.,2.)
    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    X, Y = m([min_cs_lon, max_cs_lon, max_cs_lon, min_cs_lon,min_cs_lon], [max_cs_lat,max_cs_lat,min_cs_lat, min_cs_lat,max_cs_lat])
    m.plot(X,Y, color='k', markersize=300,zorder=50)

    ax2 = plt.subplot2grid((2,3), (0,1), colspan=1)
    zoom_rescale = abs(0.6*(max(profile.length)/111100))
    m = Basemap(projection='merc', llcrnrlat=med_cs_lat-zoom_rescale,
                urcrnrlat=med_cs_lat+zoom_rescale,llcrnrlon=med_cs_lon-2*zoom_rescale,
                urcrnrlon=med_cs_lon+2*zoom_rescale, resolution='l') #lat_ts=-25
    m.drawmapboundary(fill_color='#99ffff')
    m.readshapefile('./shapefiles/MI_outline', 'MI_outline')
    m.readshapefile('./shapefiles/WI_outline', 'WI_outline')
    m.readshapefile('./shapefiles/MN_outline', 'MN_outline')
    patches   = []
    for info, shape in zip(m.MI_outline_info, m.MI_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax2.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))
    patches   = []
    for info, shape in zip(m.WI_outline_info, m.WI_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax2.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))
    patches   = []
    for info, shape in zip(m.MN_outline_info, m.MN_outline):
        patches.append( Polygon(np.array(shape), True) )
    ax2.add_collection(PatchCollection(patches, facecolor= '#cc9966', edgecolor='none', linewidths=1., zorder=2))

    # m.drawrivers(color='#99ffff')
    # m.drawcoastlines()
    # m.drawmapboundary(fill_color='#99ffff')
    # m.fillcontinents(color='#cc9966',lake_color='#99ffff')
    # parallels = np.arange(-90,90,2.)
    # m.drawparallels(parallels,labels=[1,0,0,0],fontsize=10)
    # meridians = np.arange(0.,360.,2.)
    # m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
    X, Y = m([profile['lon'][min_ind], profile['lon'][max_ind]],
             [profile['lat'][min_ind], profile['lat'][max_ind]])
    m.plot(X, Y, linestyle='dotted', linewidth=2, zorder=50)
    if locations:
        for location in locations:
            X, Y = m(loc_coor[location][1], loc_coor[location][0])
            m.scatter(X,Y, label=location, color='r', zorder=52)
    X, Y = m(profile['lon'][min_ind], profile['lat'][min_ind])
    plt.text(X, Y, 'A', bbox=dict(facecolor='white', edgecolor='white'), horizontalalignment='center', zorder=51)
    X, Y = m(profile['lon'][max_ind], profile['lat'][max_ind])
    plt.text(X, Y, 'A\'', bbox=dict(facecolor='white', edgecolor='white'), horizontalalignment='center', zorder=51)
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
    plt.xlim(min(profile.length), max(profile.length))
    plt.tight_layout()

    ax4 = plt.subplot2grid((2,3), (1, 0), colspan=3)
    plt.plot(profile.length.tolist(), profile.elevation.tolist())
    contact_strat = {}
    if reset_strat0 is not None:
        strat = -reset_strat0
    else:
        strat=0
    structures = structures.loc[structures['length']>=min_length].loc[structures['length']<=max_length]
    structures.sort_values(by='length', inplace=True)
    structures.reset_index(inplace=True, drop=True)
    for n in range(len(structures)):
        x_linspace = np.linspace(structures['length'][n], structures['length'][n]+(structures['elevation'][n]-ymin)/np.tan(np.deg2rad(abs(structures['corr_dip'][n]))))
        if n==0:
            try:
                dip_panel = np.array([structures['dip'][n], structures['dip'][n+1]])
            except:
                dip_panel = np.array([structures['dip'][n]])
        else:
            try:
                dip_panel = np.array([structures['dip'][n], structures['dip'][n-1]])
            except:
                dip_panel = np.array([structures['dip'][n]])
        avg_dip = np.average(dip_panel)
        y = structures['slope'][n]*(x_linspace-structures['length'][n])+structures['elevation'][n]
        if n>0:
            sep=structures['length'][n]-structures['length'][n-1]
        else:
            sep=0
        plt.plot(x_linspace, y, 'r')
        if id_measurements==True:
            site_name = str(structures['id'][n])
            plt.text(x_linspace[0], y[0], site_name, color='black', weight='semibold', bbox=dict(facecolor='white', edgecolor='white', alpha=0.8), horizontalalignment='center')
        dip_d_angle = calc_dip_trend_diff(profile_trend,structures['dip_d'][n])
        prior_strat = strat
        if n==0:
            elevation_change=0
        else:
            elevation_change = float(structures['elevation'][n]-structures['elevation'][n-1])
        strat += calc_strat(sep, dip_d_angle, avg_dip, elevation_change)
        extrap_strat = strat
        #print(sep, dip_d_angle, avg_dip, elevation_change)
        # strat += (sep/np.cos(np.deg2rad(dip_d_angle))*np.sin(np.deg2rad(abs(avg_dip))))+(elevation_change*np.cos(np.deg2rad(abs(avg_dip))))
        lab_height = 0.1*(ylim-max(y)) + 0.5*(ylim-max(y))*(np.random.rand(1)+0.1*(n%6))#(30 + (n * 10)%60)*
        plt.text(structures['length'][n], max(y)+lab_height, '{:}'.format(int(strat)), bbox=dict(facecolor='white', edgecolor='white'), horizontalalignment='center')
        plt.vlines(structures['length'][n], structures['elevation'][n], max(y)+lab_height-4, linestyles='dotted')
        #record strat heights in all_strat_heights dictionary
        try:
            all_strat_heights[str(structures['id'][n])] = int(strat)
        except:
            pass
        # record the strat height of contacts (if entered)
        if contacts:
            for contact in contacts:
                current_elevation = profile.loc[profile['length']==find_nearest(profile['length'], contacts[contact])]['elevation'].tolist()[0]
                prior_elevation = structures['elevation'][0]
                if n==0 and contacts[contact]<=structures['length'][n]:
                    elevation_change = profile.loc[profile['length']==find_nearest(profile['length'], contacts[contact])]['elevation'].tolist()[0]-structures['elevation'][0]
                    contact_length_diff = contacts[contact]-structures['length'][0]
                    contact_strat[contact] = prior_strat+calc_strat(contact_length_diff, dip_d_angle, avg_dip, elevation_change)
                #     if contacts[contact]<=structures['length'][n] and contacts[contact]>=profile['length'][0]:
                #         contact_strat[contact] = extrap_strat+(contacts[contact]-structures['length'][n])*np.sin(np.deg2rad(abs(structures['corr_dip'][n])))
                elif contacts[contact]<=structures['length'][n] and contacts[contact]>=structures['length'][n-1]:
                    elevation_change = profile.loc[profile['length']==find_nearest(profile['length'], contacts[contact])]['elevation']-structures['elevation'][n-1]
                    contact_length_diff = contacts[contact]-structures['length'][n-1]
                    contact_strat[contact] = prior_strat+calc_strat(contact_length_diff, dip_d_angle, avg_dip, elevation_change)
        if n>0 and n<len(structures)-1:
            if structures['length'][n]-structures['length'][n-1]>extrap_range:
                extrap_strat=prior_strat
                for extrap in range(int(structures['length'][n-1])+extrap_range, int(structures['length'][n]), extrap_range):
                    current_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'].values
                    prior_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap-extrap_range)]['elevation'].values
                    elevation_change = current_elevation - prior_elevation
                    extrap_strat += calc_strat(extrap_range, dip_d_angle, avg_dip, elevation_change)
                    # extrap_strat += extrap_range/np.cos(np.deg2rad(dip_d_angle))*np.sin(np.deg2rad(avg_dip))
                    # text_height = 30 + (extrap * 10)%60
                    plt.text(extrap,
                             profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height, '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                    plt.vlines(extrap,
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height-5,
                               linestyles='dotted')
        elif n==len(structures)-1:
            if max_length-structures['length'][n]>extrap_range:
                # extrap_strat = strat
                for extrap in range(int(structures['length'][n])+extrap_range, max_length, extrap_range):
                    current_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'].values
                    prior_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap-extrap_range)]['elevation'].values
                    elevation_change = current_elevation - prior_elevation
                    extrap_strat += calc_strat(extrap_range, dip_d_angle, avg_dip, elevation_change)
                    # extrap_strat += extrap_range/np.cos(np.deg2rad(dip_d_angle))*np.sin(np.deg2rad(avg_dip))
                    plt.text(extrap,
                             profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height,
                             '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                    plt.vlines(extrap,
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height-5,
                               linestyles='dotted')
            # calculate total strat of profile
            farthest_extrap=(max_length-int(structures['length'][n])+extrap_range)//extrap_range
            current_elevation = profile.loc[profile['length']==find_nearest(profile['length'], max(profile.length))]['elevation'].values
            prior_elevation = profile.loc[profile['length']==find_nearest(profile['length'], farthest_extrap)]['elevation'].values
            elevation_change = current_elevation - prior_elevation
            total_strat = extrap_strat + calc_strat(max(profile.length)-farthest_extrap, dip_d_angle, avg_dip, elevation_change)
        if n==0:
            extrap_strat=0
            if (structures['length'][n] - min_length)>extrap_range:
                for extrap in range(int(min(structures.length))-extrap_range, min_length, -extrap_range):
                    current_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'].values
                    prior_elevation = profile.loc[profile['length']==find_nearest(profile['length'], extrap+extrap_range)]['elevation'].values
                    elevation_change = current_elevation - prior_elevation
                    # print(extrap_strat)
                    extrap_strat += calc_strat(-extrap_range, dip_d_angle, avg_dip, elevation_change)
                    # print(-extrap_range, dip_d_angle, avg_dip, elevation_change)
                    # print(extrap_strat)
                    # extrap_strat += -extrap_range*np.cos(np.deg2rad(dip_d_angle))*np.sin(np.deg2rad(abs(structures['dip'][0])))
                    plt.text(extrap,
                             profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height,
                             '{0}'.format(int(extrap_strat)), horizontalalignment='center')
                    plt.vlines(extrap,
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation'],
                               profile.loc[profile['length']==find_nearest(profile['length'], extrap)]['elevation']+lab_height-5,
                               linestyles='dotted')
    if contacts:
        for contact in contacts:
            xline = float(contacts[contact])
            contact_prof_index = profile.loc[profile['length']==find_nearest(profile['length'], xline)].index.tolist()[0]
            contact_str_index = structures.loc[structures['length']==find_nearest(structures['length'], xline)].index.tolist()[0]
            if xline <= structures['length'][contact_str_index]:
                try:
                    dip_panel = np.array([structures['corr_dip'][contact_str_index], structures['corr_dip'][contact_str_index-1]])
                except:
                    dip_panel = np.array([structures['corr_dip'][contact_str_index]])
            elif xline >= structures['length'][contact_str_index]:
                try:
                    dip_panel = np.array([structures['corr_dip'][contact_str_index], structures['corr_dip'][contact_str_index+1]])
                except:
                    dip_panel = np.array([structures['corr_dip'][contact_str_index]])
            avg_dip = abs(np.average(dip_panel))
            x_contact = np.linspace(xline, xline+(profile['elevation'][contact_prof_index]-ymin)/np.tan(np.deg2rad(avg_dip)))
            y_contact = -abs(np.tan(np.deg2rad(avg_dip)))*(x_contact-profile['length'][contact_prof_index])+profile['elevation'][contact_prof_index]
            ymax = profile['elevation'][contact_prof_index]
            contact_color = np.array([0.,0.4,0.4])
            # contact_color = np.random.rand(3)
            plt.plot(x_contact,y_contact,color=contact_color, linewidth=2.0, linestyle='dashed', label=contact, alpha=0.8)
            y_place = np.random.randint(ymax+0.1*(ylim-ymin),ylim-0.1*(ylim-ymin))
            plt.text(xline, y_place+5, '{:d}'.format(int(contact_strat[contact])), color='w',weight='semibold', bbox=dict(facecolor=contact_color, edgecolor='none'),horizontalalignment='center')
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

    ax4.text(0, 1, 'A',
             color='k', weight='semibold', size='large', transform=ax4.transAxes,
             bbox=dict(facecolor='white', edgecolor='white'), verticalalignment='top',zorder=50)
    ax4.text(1, 1, 'A\'',
             color='k', weight='semibold', size='large',  transform=ax4.transAxes,
             bbox=dict(facecolor='white', edgecolor='white'), verticalalignment='top', zorder=50)
    # ax4.text(0, 1, '{:.5}\n{:.5}'.format(profile['lat'][min_ind],
    #                                      profile['lon'][min_ind]),
    #          color='r', transform=ax4.transAxes,
    #          bbox=dict(facecolor='white', edgecolor='white'), verticalalignment='top',zorder=50)
    # ax4.text(1, 1, '{:.5}\n{:.5}'.format(profile['lat'][max_ind],
    #                                      profile['lon'][max_ind]),
    #          color='r', transform=ax4.transAxes,
    #          bbox=dict(facecolor='white', edgecolor='white'), verticalalignment='top', zorder=50)

    plt.xlabel('meters')
    plt.ylabel('meters')

    plt.xlim(min_length,max_length)
    plt.ylim(ymin,ylim)

    plt.title('Bearing = {:.4}°, Total Strat = {:.4} km'.format(calc_cs_trend(profile), int(total_strat)/1000))

    plt.legend()

    if save==True:
        plt.savefig(str(save_name))

    plt.show()

    if return_all_strat_heights:
        return all_strat_heights


def make_cross_section(profile_data, str_data, return_structures=False, **kwargs):
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
    if return_structures:
        return cs_figure(profile, structures, **kwargs), structures
    return cs_figure(profile, structures, **kwargs)
