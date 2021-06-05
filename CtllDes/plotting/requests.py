import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_comm(comm_df, mask = None, deg=False, total_time=False):
    """Plot Communicatoin data.
    
    Parameters
    ----------
    comm_df : pandas.DataFrame
        Output of CtllDes.requests.comm.default_comm_data(...) 
    mask : float
        Minimum angle required to establish communication
    deg : boolean, optional
        set to degrees, if mask is used it must match the unit used.
    total_time : bollean, optional
        report in figure legend total time over ground station.
    """
    
    if not deg:
        elevation = comm_df['elevation'].to_numpy(dtype=np.float64)
    else:
        elevation = comm_df['elevation'].to_numpy(dtype=np.float64)*180/np.pi
    
    time = comm_df['time']/3600
    
    fig, ax = plt.subplots(figsize=(10,10))
    ax.plot(time, elevation, linewidth=1, color='black', label='Angle to ground station')
    ax.grid()
    if mask != None:
        elevation_mask = [e for e in elevation if e > mask]
        time_over = len(elevation_mask)*(time[1]-time[0])*3600 #very nasty

        if not deg:
            ax.hlines(mask, xmin=0, xmax=max(time), label=f'Mask Angle = {mask}, time over mask = {time_over/60:.2f} minutes')
        else:
            ax.hlines(mask, xmin=0, xmax=max(time),label=f'Mask Angle = {mask}, time over mask = {time_over/60:.2f} minutes ')
            
    ax.set_xlabel("Tiempo [h]")
    if not deg:
        ax.set_ylabel("Elevation [rad]")
    else:
        ax.set_ylabel("Elevation [°]")        
    
    ax.legend(loc=1)


def plot_heat_map(covdf,
                  tgts,
                  target_name = 'Not specified',
                  accumulated=True,
                  response=True,
                  average_time_gap=True,
                  revisit=True):

    """Plot heat map of coverage data.

    Parameters
    ----------
    covdf : ~pandas.DataFrame
        Coverage dataframe, output of CtllDes.request.coverage.Coverages
    tgts : ~CtllDes.targets.targets.Targets
        Targets of analysis
    target_name : string, optional
        description of Target Name, will be setted in the title of the figure.
    accumulated : boolean, optional
        Required to plot Total accumulated time of coverage
    response : boolean, optional
        Required to plot Mean Response time
    average_time_gap : boolean, optional
        Required to plot Average Time Gap
    revisit : boolean, optional
        Required  to plot mean revisit time
        
    """

    target_lons = np.round(tgts.lons, decimals=2)
    target_lats = np.round(tgts.lats, decimals=2)

    roll_accum = covdf['accumulated'].to_numpy(dtype=float)/3600
    #roll_accum /= max(roll_accum)

    response_time = covdf['response time'].to_numpy()/3600
    response_time = 1/response_time
    response_time -= min(response_time)
    response_time /= max(response_time)
        
    roll_avg = covdf['average time gap'].to_numpy(dtype=float)/3600
    roll_avg = 1/roll_avg
    roll_avg -= min(roll_avg)
    roll_avg /= max(roll_avg) 

    revisit_time = covdf['mean gap dark'].to_numpy(dtype=np.float64)/3600
        
    if accumulated:
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes()
        
        df = pd.DataFrame.from_dict(np.array([target_lons,target_lats,roll_accum]).T)
        df.columns = ['Longitude [°]','Latitude [°]','Accumulated time of view [hours]']
        df['Accumulated time of view [hours]'] = pd.to_numeric(df['Accumulated time of view [hours]'])
        pivotted= df.pivot('Latitude [°]','Longitude [°]','Accumulated time of view [hours]')
        sns.heatmap(pivotted, cmap='magma', ax=ax)
        ax.set_title(f"Accumulated time of view over {target_name}")

    if response:
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes()
        
        df = pd.DataFrame.from_dict(np.array([target_lons,target_lats,response_time]).T)
        df.columns = ['Longitude [°]','Latitude [°]','Mean Response Time [1/hours]']
        df['Mean Response Time [1/hours]'] = pd.to_numeric(df['Mean Response Time [1/hours]'])
        pivotted= df.pivot('Latitude [°]','Longitude [°]','Mean Response Time [1/hours]')
        sns.heatmap(pivotted, cmap='magma', ax=ax)
        ax.set_title(f"Mean Response Time [1/hours] over {target_name}")

    if revisit:
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes()
        
        df = pd.DataFrame.from_dict(np.array([target_lons,target_lats,revisit_time]).T)
        df.columns = ['Longitude [°]','Latitude [°]','Mean Revisit Time [hours]']
        df['Mean Revisit Time [hours]'] = pd.to_numeric(df['Mean Revisit Time [hours]'])
        pivotted= df.pivot('Latitude [°]','Longitude [°]','Mean Revisit Time [hours]')
        sns.heatmap(pivotted, cmap='magma', ax=ax)
        ax.set_title(f"Mean Revisit Time [hours] over {target_name}")        

    if average_time_gap:
        fig = plt.figure(figsize=(10,10))
        ax = plt.axes()
        
        df = pd.DataFrame.from_dict(np.array([target_lons,target_lats,roll_accum]).T)
        df.columns = ['Longitude [°]','Latitude [°]','Average Time Gap [1/hours]']
        df['Average Time Gap [1/hours]'] = pd.to_numeric(df['Average Time Gap [1/hours]'])
        pivotted= df.pivot('Latitude [°]','Longitude [°]','Average Time Gap [1/hours]')
        sns.heatmap(pivotted, cmap='magma', ax=ax)
        ax.set_title(f"Average Time Gap [1/hours] over {target_name}")


def plot_Targets_merit(rollcovdf,
                       tgts,
                       target_name = 'Not specified',
                       accumulated=True,
                       response=True,
                       average_time_gap=True,
                       revisit=True,
                       use_3d=False):
    """Plot heat map of coverage data.

    Parameters
    ----------
    covdf : ~pandas.DataFrame
        Coverage dataframe, output of CtllDes.request.coverage.Coverages
    tgts : ~CtllDes.targets.targets.Targets
        Targets of analysis
    target_name : string, optional
        description of Target Name, will be setted in the title of the figure.
    accumulated : boolean, optional
        Required to plot Total accumulated time of coverage
    response : boolean, optional
        Required to plot Mean Response time
    average_time_gap : boolean, optional
        Required to plot Average Time Gap
    revisit : boolean, optional
        Required  to plot mean revisit time
    use_3d : boolean, optional
        Plot trisfurf in 3d.    
    """


    target_lons = tgts.lons
    target_lats = tgts.lats

    roll_accum = rollcovdf['accumulated'].to_numpy(dtype=float)/3600
    #roll_accum /= max(roll_accum)

    response_time = rollcovdf['response time'].to_numpy()/3600
    response_time = 1/response_time
    response_time -= min(response_time)
    response_time /= max(response_time)
        
    roll_avg = rollcovdf['average time gap'].to_numpy(dtype=float)/3600
    roll_avg = 1/roll_avg
    roll_avg -= min(roll_avg)
    roll_avg /= max(roll_avg) 
    
    revisit_time = rollcovdf['mean gap dark'].to_numpy(dtype=np.float64)/3600
        
    if accumulated:
        fig = plt.figure(figsize=(10,10))
        
        if use_3d:
            ax.plot_trisurf(target_lons, target_lats, roll_accum,
                       antialiased=False, cmap='afmhot')
            ax.scatter(target_lons,target_lats, np.zeros(len(target_lons)),
                   s=100,c='k')
        else:
            ax = plt.axes()
            tcs = ax.tricontourf(target_lons, target_lats, roll_accum)
            cb = plt.colorbar(tcs)
            cb.ax.set_ylabel("Accumulated time of view [hours]")

        
        ax.set_title(f"Accumulated time of view for {target_name} ")
        
        ax.set_xlabel("longitude [°]")
        ax.set_ylabel("latitude [°]")
        if use_3d:
            ax = fig.add_subplot(projection='3d')
            ax.set_zlabel("Accumulated time of view [hours]")
            
        ax.set_xlim(min(target_lons),max(target_lons))
        ax.set_ylim(min(target_lats),max(target_lats))


        
    if response:
        fig = plt.figure(figsize=(10,10))
        
        if use_3d:
            ax = fig.add_subplot(projection='3d')
            ax.plot_trisurf(target_lons, target_lats, response_time,
                       antialiased=False, cmap='afmhot')
            ax.scatter(target_lons,target_lats, np.zeros(len(target_lons)),
                   s=100,c='k')
        else:
            ax = plt.axes()
            tcs = ax.tricontourf(target_lons, target_lats, response_time)
            cb = plt.colorbar(tcs)
            cb.ax.set_ylabel("Respose time [1/hours]")
            
        ax.set_title(f"Response Time for {target_name} ")
        
        ax.set_xlabel("longitude [°]")
        ax.set_ylabel("latitude [°]")
        if use_3d:
            ax.set_zlabel("Respose time [1/hours]")
            
        ax.set_xlim(min(target_lons),max(target_lons))
        ax.set_ylim(min(target_lats),max(target_lats))
        

        
    if average_time_gap:        
        fig = plt.figure(figsize=(10,10))
        
        if use_3d:
            ax = fig.add_subplot(projection='3d')
            ax.plot_trisurf(target_lons, target_lats, roll_avg,
                       antialiased=False, cmap='afmhot')
            ax.scatter(target_lons,target_lats, np.zeros(len(target_lons)),
                   s=100,c='k')
        else:
            ax = plt.axes()
            tcs = ax.tricontourf(target_lons, target_lats, roll_avg)
            cb = plt.colorbar(tcs)
            cb.ax.set_ylabel("Average time Gap for [1/hours]")

        ax.set_title(f"Average time Gap for {target_name} ")
        
        ax.set_xlabel("longitude [°]")
        ax.set_ylabel("latitude [°]")
        if use_3d:
            ax.set_zlabel("Average time Gap for [1/hours]")
            
        ax.set_xlim(min(target_lons),max(target_lons))
        ax.set_ylim(min(target_lats),max(target_lats))
        

       
    if revisit:
        fig = plt.figure(figsize=(10,10))
        
        if use_3d:
            ax = fig.add_subplot(projection='3d')
            ax.plot_trisurf(target_lons, target_lats, revisit_time,
                       antialiased=False, cmap='afmhot')
            ax.scatter(target_lons,target_lats, np.zeros(len(target_lons)),
                   s=100,c='k')
        else:
            ax = plt.axes()
            tcs = ax.tricontourf(target_lons, target_lats, revisit_time)
            cb = plt.colorbar(tcs)
            cb.ax.set_ylabel("Mean Revisit Time [hours]")
        
        ax.set_title(f"Mean Revisit Time {target_name} ")
        
        ax.set_xlabel("longitude [°]")
        ax.set_ylabel("latitude [°]")
        if use_3d:
            ax.set_zlabel("Mean Revisit Time [hours]")
            
        ax.set_xlim(min(target_lons),max(target_lons))
        ax.set_ylim(min(target_lats),max(target_lats))

