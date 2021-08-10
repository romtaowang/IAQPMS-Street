#!/usr/bin/python3
# -*- coding: utf-8 -*-

import io
import os,sys
import math,re
import numpy as np
import datetime,time
import scipy.io as scio
from scipy.sparse import coo_matrix
from mpi4py import MPI

import ef_equation_china
import sys

reload(sys)
sys.setdefaultencoding('utf8')

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size= comm.Get_size()
if __name__ == "__main__":
    work_path = '../../'
    base_info_path= work_path + 'base_info/'
    road_map_file= base_info_path + 'post_map.csv'
    road_inring_file= base_info_path + 'id.csv'
    data_path= work_path +'/calculation/'
    output_path= work_path +'/calculation/emission_hour/test_mpi6/'
    
    ### input data can be updated
    
    vehicle_fleet_flag_segment=['Small','Medium','Large','Small','Medium','Large','Motorcycles 4-stroke 250 - 750 cm³','Urban Buses Standard 15 - 18 t','Small']
    
    
    vehicle_fleet_flag=['LDV','MDV','HDV','LDT','MDT','HDT','MC','Bus','Taxi']
    vehicle_fleet=np.loadtxt(base_info_path + 'vehicle_fleet.csv',delimiter = ',', skiprows = 1,usecols=(1,))
    vehicle_fleet2=np.loadtxt(base_info_path + 'vehicle_percentage_Beijing_out.csv',delimiter = ',', skiprows = 1,usecols=(1,))
    print('vehicle read ok') 
    fuel_type_flag=['Petrol','Diesel','LPG']
    fuel_type_percentage=np.loadtxt(base_info_path + 'fuel_type_percentage.csv',delimiter = ',', skiprows = 1,usecols=(1,2,3))
    fuel_type_percentage=fuel_type_percentage/100
    
    emission_standard_flag=['Pre I','CHINA I','CHINA II','CHINA III','CHINA IV','CHINA V']
    
    
    emission_standard_percentage=np.loadtxt(base_info_path + 'emission_standard_percentage.csv',delimiter = ',', skiprows = 1,usecols=(1,2,3,4,5,6))
    emission_standard_percentage=emission_standard_percentage/100
    ###
    
    
    species_flag=['CO','HC','NOx','PM25','PM10']
    
    
    f_basemap = io.open(road_map_file, 'r', encoding='utf-8')
    
    x_s=[];y_s=[];x_e=[];y_e=[];base_coord=[];
    fclass=[];oneway=[];street_length=[];
     
    for line in f_basemap.readlines()[1:]:
        # print(line)
        x_s.append(line.split(',')[5])
        y_s.append(line.split(',')[6])
        x_e.append(line.split(',')[7])
        y_e.append(line.split(',')[8])
        base_coord.append(re.split(',|\n',line)[5:9])
        fclass.append(line.split(',')[9])
        oneway.append(line.split(',')[10])
        street_length.append(re.split(',|\n',line)[26])
    f_basemap.close()
    print('basemap ok')
    in_basemap = io.open(road_inring_file, 'r', encoding='utf-8')
    in_s=[]
    for line2 in in_basemap.readlines()[1:]:
        in_s.append(line2.split(',')[1])
    in_basemap.close()
    ttt=np.array(in_s)
    ttt = ttt.astype(np.int)
    print('inbasemap ok')
    print(ttt[10])
    f_emission_factor=io.open(base_info_path + 'On Road Hot Emission Factor Chinese.csv', 'r', encoding='utf-8')
    
    emission_factor=[];emission_factor_flag=[]
    emission_factor_flag.append(re.split(',|\n',f_emission_factor.readline()))
    
    for line in f_emission_factor.readlines():
        emission_factor.append(re.split(',|\n',line))
    
    f_emission_factor.close()
    
    #############
    ### for test
    #l=n=h=s=vf=ft=es=ef=0
    #############
    
    high_altitude = False  # whether the altitude is higher than 1500m
    
    ### for evporation  in hour
    ef1 = 11.6                # evporation emission factor during driving
    ef2 = 6.5 / 24 * (3/60)   # evporation emission factor during parking
    ###
    
    
    start_day=datetime.datetime(2021,7,1)
    end_day=datetime.datetime(2021,7,31)
    
    span_day=(end_day-start_day).days+1
    
    
    
    
    for n in range(span_day):
        start_time_cal = time.time()
        cal_day = start_day + datetime.timedelta(days=n)
        print(cal_day)
        flow_data=np.loadtxt(data_path + '/flow_hour/post_flow/' \
                              +str(cal_day.year).zfill(4)+'-'+str(cal_day.month).zfill(2)+'-'+str(cal_day.day).zfill(2)+'_hour_flow.txt', \
                              delimiter = ',', skiprows = (1))
    
        speed_data = np.loadtxt(data_path + '/speed_hour/post_speed/' \
                                + str(cal_day.year).zfill(4) + '-' + str(cal_day.month).zfill(2) + '-' + str(cal_day.day).zfill(2) + '_hour_speed.txt', \
                                delimiter=',', skiprows=(1))
        #if rank == 0:
        total_emission_output=np.zeros([int(math.ceil(float(len(base_coord))/size))*size, 24, len(species_flag)],'float')
        #total_emission = np.zeros([int(math.ceil(float(len(base_coord))/size))*size, 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        #total_emission=coo_matrix((int(math.ceil(float(len(base_coord))/size))*size, 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)),dtype=np.float).toarray()
        print('total_emission size 0 is')
        print(int(math.ceil(float(len(base_coord))/size))*size)
        #hef = np.zeros([int(math.ceil(float(len(base_coord))/size))*size, 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        street_list=len(base_coord)
        h_list=24
        species_list=len(species_flag)
        vehicle_fleet_list=len(vehicle_fleet_flag)
        fuel_list=len(fuel_type_flag)
        emission_list=len(emission_standard_flag)
        #else:
        #    street_list=0
        #    h_list=0
        #    species_list=0
        #    vehicle_fleet_list=0
        #    fuel_list=0
        #    emission_list=0
        #    total_emission = comm.recv(source=0)
        #    total_emission_output = comm.recv(source=0)
    
        #street_list=comm.bcast(street_list,root=0)
        #h_list=comm.bcast(h_list,root=0)
        #species_list=comm.bcast(species_list,root=0)
        #vehicle_fleet_list=comm.bcast(vehicle_fleet_list,root=0)
        #fuel_list=comm.bcast(fuel_list,root=0)
        #emission_list=comm.bcast(emission_list,root=0)
        #total_emission = np.zeros([street_list,h_list,species_list,vehicle_fleet_list,fuel_list,emission_list],'float')
           #comm.Bcast(total_emission,root=0)
        ##total_emission=comm.bcast(total_emission if rank == 0 else None, root=0)
        num_samples=int(math.ceil(float(len(base_coord))/size))*size     #total_emission.shape[0]
        local_data_offset=np.linspace(0, num_samples, size + 1).astype('int') 
        #local_data=total_emission[local_data_offset[rank] :local_data_offset[rank + 1]]
        #total_emission_i=np.zeros([local_data.shape[0],local_data.shape[1],local_data.shape[2],local_data.shape[3],local_data.shape[4],local_data.shape[5]],'float')
        #hef_i=np.zeros([local_data.shape[0],local_data.shape[1],local_data.shape[2],local_data.shape[3],local_data.shape[4],local_data.shape[5]],'float')
        total_emission_i=np.zeros([int(math.ceil(float(len(base_coord))/size)), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        hef_i=np.zeros([int(math.ceil(float(len(base_coord))/size)), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        #if rank+1==size:
        #    total_emission_i = np.zeros([len(base_coord)-int(math.floor(len(base_coord)//size))*(size-1), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        #    hef_i = np.zeros([len(base_coord)-int(math.floor(len(base_coord)//size))*(size-1), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        #else:
        #    a=int(math.floor(len(base_coord)//size))
        #    print(str(a))
        #    total_emission_i = np.zeros([a, 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        #    hef_i = np.zeros([a, 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)],'float')
        for h in range(24):
            print('nday=' + str(n + 1) + '; hour=' + str(h))
            vehicle_fleet_number = {}
            #vehicle_fleet_percentage = {}
            vehicle_fleet_number2 = {}
            vehicle_fleet_percentage2 = {}
            vehicle_fleet_percentageback = {}
            for vff in range(len(vehicle_fleet_flag)):
                vehicle_fleet_number[vehicle_fleet_flag[vff]] = vehicle_fleet[vff]
                vehicle_fleet_number2[vehicle_fleet_flag[vff]] = vehicle_fleet2[vff]
            ### traffic control
            # for Motorcycles
            vehicle_fleet_number['MC'] = 0
            # for trunk
            if 7<=h<=22:
                vehicle_fleet_number['MDT'] = 0
                vehicle_fleet_number['HDT'] = 0
            if 7<=h<=9 or 17<=h<=20:
                vehicle_fleet_number['LDT'] = 0
            # for bus
            if 2 <= h <= 5:
                vehicle_fleet_number['Bus'] = 0
            
            ###
            total_vehicle_fleet_number = sum(i for i in vehicle_fleet_number.values())
            total_vehicle_fleet_number2 = sum(i for i in vehicle_fleet_number2.values())
            for vff in range(len(vehicle_fleet_flag)):
                vehicle_fleet_percentageback[vehicle_fleet_flag[vff]] = vehicle_fleet_number[vehicle_fleet_flag[vff]] / total_vehicle_fleet_number
                vehicle_fleet_percentage2[vehicle_fleet_flag[vff]] = vehicle_fleet_number2[vehicle_fleet_flag[vff]] / total_vehicle_fleet_number2
            #if rank+1==size:
            #    total_emission_i = np.zeros((len(base_coord)-math.floor(len(base_coord)//size)*(size-1), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)))
            #    hef_i = np.zeros((len(base_coord)-math.floor(len(base_coord)//size)*(size-1), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)))
            #else:	
            #    total_emission_i = np.zeros((math.floor(len(base_coord)//size), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)))
            #    hef_i = np.zeros((math.floor(len(base_coord)//size), 24, len(species_flag), len(vehicle_fleet_flag), len(fuel_type_flag),len(emission_standard_flag)))
            ### calculate the emission factor
            #if rank+1==size:
            #    begin_l=rank*int(math.floor(len(base_coord)//size))
            #    end_l=len(base_coord)
            #else:
            #    begin_l=rank*int(math.floor(len(base_coord)//size))
            #    end_l=(rank+1)*int(math.floor(len(base_coord)//size))
            
            begin_l=local_data_offset[rank]
            end_l=local_data_offset[rank + 1]
            #for l in range(rank*len(base_coord)//size,(rank+1)*len(base_coord)//size):
            for l in range(begin_l,end_l):
              if (l < len(base_coord)):
                # for l in range(0,11):
                #     print(l)
                vehicle_fleet_percentage = {}
                if (ttt[l]==1):
                   for vff in range(len(vehicle_fleet_flag)):
                       vehicle_fleet_percentage[vehicle_fleet_flag[vff]] = vehicle_fleet_percentageback[vehicle_fleet_flag[vff]]
                else:
                   for vff in range(len(vehicle_fleet_flag)):
                       vehicle_fleet_percentage[vehicle_fleet_flag[vff]] = vehicle_fleet_percentage2[vehicle_fleet_flag[vff]]
                
                for s in range(len(species_flag)):
                    ### for Passenger Car
                    for vf in range(len(vehicle_fleet_flag[0:3])):
                        for ft in range(len(fuel_type_flag)):
                            for es in range(len(emission_standard_flag)):
                                ### for evaporation
                                if species_flag[s] == 'HC' and fuel_type_flag[ft] == 'Petrol':
                                    total_emission_i[l-begin_l, h, s, vf, ft, es] = total_emission_i[l-begin_l, h, s, vf, ft, es] + (ef1 * float(street_length[l]) /1000 / speed_data[l,h] + ef2) * flow_data[l,h] \
                                                                           * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                                          fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es]

                                for ef in range(len(emission_factor)):
                                    ### for exhaust
                                    if emission_factor[ef][0] == 'Passenger Cars' and emission_factor[ef][1] == fuel_type_flag[ft] and \
                                            emission_factor[ef][2] == vehicle_fleet_flag_segment[vf] and emission_factor[ef][3] == emission_standard_flag[es] and \
                                            emission_factor[ef][5] == species_flag[s]:
                                        try:
                                            c_speed = ef_equation_china.speed_correction(speed_data[l,h], species_flag[s], emission_factor[ef][1], emission_standard_flag[es])

                                        except UnboundLocalError:
                                            c_speed = 1

                                        try:
                                            c_ta = ef_equation_china.ta_correction(cal_day.month, species_flag[s], emission_factor[ef][1], vehicle_fleet_flag[vf])

                                        except UnboundLocalError:
                                            c_ta = 1

                                        try:
                                            c_rh = ef_equation_china.rh_correction(cal_day.month, species_flag[s], emission_factor[ef][1])

                                        except UnboundLocalError:
                                            c_rh = 1

                                        if high_altitude:
                                            try:
                                                c_height = ef_equation_china.height_correction(species_flag[s],emission_factor[ef][1],vehicle_fleet_flag[vf])

                                            except UnboundLocalError:
                                                c_height = 1

                                        else:
                                            c_height = 1

                                        # c_speed = c_ta = c_rh = c_height = 1
                                        hef_i[l-begin_l,h,s,vf,ft,es] = ef_equation_china.ef_equation_china(float(emission_factor[ef][6]), c_speed, c_ta, c_rh, c_height)
                                        #print(vehicle_fleet_flag[vf], speed_data[l,h], species_flag[s], fuel_type_flag[ft],emission_standard_flag[es], c_speed)
                                        #add by Tao Wang, ref: B. Y. Jing,ACP,for Bejing
                                        #if species_flag[s] == 'NOx':
                                        #   hef_back=hef[l,h,s,vf,ft,es]
                                        #   if vf == 0:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*0.7067-0.0604
                                        #   elif vf == 1:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*0.5659-0.3019
                                        #   else:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*1.6664-7.6978
                                        #   if hef[l,h,s,vf,ft,es] <= 0:
                                        #      hef[l,h,s,vf,ft,es] = hef_back
                                        #
                                        total_emission_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * \
                                                              flow_data[l,h] * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                              fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es] * \
                                                              float(street_length[l]) /1000
                                        break
                    ### for Truck and motorcycle
                    for vf in range(len(vehicle_fleet_flag[0:3]), len(vehicle_fleet_flag[0:7])):
                        emission_standard_flag_truck = emission_standard_flag
                        if vehicle_fleet_flag[vf] == 'MC':
                            vehicle_fleet_flag_truck = 'L-Category'
                        else:
                            vehicle_fleet_flag_truck = 'Heavy Duty Trucks'
    
                        for ft in range(len(fuel_type_flag)):
                            for es in range(len(emission_standard_flag_truck)):
                                ### for evaporation
                                if species_flag[s] == 'HC' and fuel_type_flag[ft] == 'Petrol':
                                    total_emission_i[l-begin_l, h, s, vf, ft, es] = total_emission_i[l-begin_l, h, s, vf, ft, es] + (ef1 * float(street_length[l]) /1000 / speed_data[l,h] + ef2) * flow_data[l,h] \
                                                                           * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                                          fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es]
                                for ef in range(len(emission_factor)):
                                    ### for exhaust
                                    if emission_factor[ef][0] == vehicle_fleet_flag_truck and emission_factor[ef][1] == fuel_type_flag[ft] and \
                                            emission_factor[ef][2] == vehicle_fleet_flag_segment[vf] and emission_factor[ef][3] == emission_standard_flag_truck[es] and \
                                            emission_factor[ef][5] == species_flag[s]:
                                        try:
                                            c_speed = ef_equation_china.speed_correction(speed_data[l,h], species_flag[s], emission_factor[ef][1], emission_standard_flag[es])
                                        except UnboundLocalError:
                                            c_speed = 1
                                        try:
                                            c_ta = ef_equation_china.ta_correction(cal_day.month, species_flag[s], emission_factor[ef][1], vehicle_fleet_flag[vf])
                                        except UnboundLocalError:
                                            c_ta = 1
                                        try:
                                            c_rh = ef_equation_china.rh_correction(cal_day.month, species_flag[s], emission_factor[ef][1])
                                        except UnboundLocalError:
                                            c_rh = 1
                                        if high_altitude:
                                            try:
                                                c_height = ef_equation_china.height_correction(species_flag[s],emission_factor[ef][1],vehicle_fleet_flag[vf])
                                            except UnboundLocalError:
                                                c_height = 1
                                        else:
                                            c_height = 1
                                        # c_speed = c_ta = c_rh = c_height = 1
                                        hef_i[l-begin_l,h,s,vf,ft,es] = ef_equation_china.ef_equation_china(float(emission_factor[ef][6]), c_speed, c_ta, c_rh, c_height)
                                        #add by Tao Wang, ref: B. Y. Jing,ACP,for Bejing
                                        #if species_flag[s] == 'NOx':
                                        #   hef_back=hef[l,h,s,vf,ft,es]
                                        #   if vf == 3:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*0.8172-1.3015
                                        #   elif vf == 4:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*1.1241-3.2323
                                        #   elif vf == 5:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]*0.8757-1.5597
                                        #   else:
                                        #      hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es]
                                        #   if hef[l,h,s,vf,ft,es] <= 0:
                                        #      hef[l,h,s,vf,ft,es] = hef_back
                                        #
                                        #print(vehicle_fleet_flag[vf], speed_data[l,h], species_flag[s], fuel_type_flag[ft],emission_standard_flag[es], c_speed)
                                        total_emission_i[l-begin_l,h,s,vf,ft,es]=hef_i[l-begin_l,h,s,vf,ft,es] * \
                                                              flow_data[l,h] * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                              fuel_type_percentage[vf, ft] * emission_standard_percentage[vf,es] * \
                                                              float(street_length[l]) /1000
                                        break
                    ### for Bus
                    for vf in range(len(vehicle_fleet_flag[0:7]),len(vehicle_fleet_flag[0:8])):
                        for ft in range(len(fuel_type_flag)):
                            if fuel_type_flag[ft] == 'LPG':
                                fuel_type_flag_bus = 'Diesel'
                            else:
                                fuel_type_flag_bus = fuel_type_flag[ft]
                            for es in range(len(emission_standard_flag)):
                                ### for evaporation
                                if species_flag[s] == 'HC' and fuel_type_flag_bus == 'Petrol':
                                    total_emission_i[l-begin_l, h, s, vf, ft, es] = total_emission_i[l-begin_l, h, s, vf, ft, es] + (ef1 * float(street_length[l]) /1000 / speed_data[l,h] + ef2) * flow_data[l,h] \
                                                                           * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                                          fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es]
                                for ef in range(len(emission_factor)):
                                    ### for exhaust
                                    if emission_factor[ef][0] == 'Buses' and emission_factor[ef][1] == fuel_type_flag_bus and \
                                            emission_factor[ef][2] == vehicle_fleet_flag_segment[vf] and emission_factor[ef][3] == emission_standard_flag[es] and \
                                            emission_factor[ef][5] == species_flag[s]:
                                        try:
                                            c_speed = ef_equation_china.speed_correction(speed_data[l,h], species_flag[s], emission_factor[ef][1], emission_standard_flag[es])
                                        except UnboundLocalError:
                                            c_speed = 1
                                        try:
                                            c_ta = ef_equation_china.ta_correction(cal_day.month, species_flag[s], emission_factor[ef][1], vehicle_fleet_flag[vf])
                                        except UnboundLocalError:
                                            c_ta = 1
                                        try:
                                            c_rh = ef_equation_china.rh_correction(cal_day.month, species_flag[s], emission_factor[ef][1])
                                        except UnboundLocalError:
                                            c_rh = 1
                                        if high_altitude:
                                            try:
                                                c_height = ef_equation_china.height_correction(species_flag[s],emission_factor[ef][1],vehicle_fleet_flag[vf])
                                            except UnboundLocalError:
                                                c_height = 1
                                        else:
                                            c_height = 1
                                        # c_speed = c_ta = c_rh = c_height = 1
                                        hef_i[l-begin_l,h,s,vf,ft,es] = ef_equation_china.ef_equation_china(float(emission_factor[ef][6]), c_speed, c_ta, c_rh, c_height)
                                        #print(vehicle_fleet_flag[vf], speed_data[l,h], species_flag[s], fuel_type_flag[ft],emission_standard_flag[es], c_speed)
                                        if species_flag[s] == 'NOx':
                                           hef_back=hef_i[l-begin_l,h,s,vf,ft,es]
                                           hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es]*8.4436-8.66
                                           if hef_i[l-begin_l,h,s,vf,ft,es] <= 0:
                                              hef_i[l-begin_l,h,s,vf,ft,es] = hef_back
                                        ### reference S. Zhang, Y. Wu, H.liu et al, AE, 2013
                                        if fuel_type_flag[ft] == 'LPG':
                                            if species_flag[s] == 'CO':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 0.11
                                            elif species_flag[s] == 'HC':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 0.64
                                            #elif species_flag[s] == 'NOx':
                                            #    hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es] * 1.71
                                            elif species_flag[s] == 'PM25' or species_flag[s] == 'PM10':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 0.1
                                            else:
                                                hef_i[l-begin_l,h,s,vf,ft,es] = 0
                                        total_emission_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * \
                                                              flow_data[l,h] * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                              fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es] * \
                                                              float(street_length[l]) /1000
                                        break
                    ### for taxi
                    for vf in range(len(vehicle_fleet_flag[0:8]),len(vehicle_fleet_flag[0:9])):
                        for ft in range(len(fuel_type_flag)):
                            if fuel_type_flag[ft] == 'LPG':
                                fuel_type_flag_taxi = 'Petrol'
                            else:
                                fuel_type_flag_taxi = fuel_type_flag[ft]
                            for es in range(len(emission_standard_flag)):
                                ### for evaporation
                                if species_flag[s] == 'HC' and fuel_type_flag_taxi == 'Petrol':
                                    total_emission_i[l-begin_l, h, s, vf, ft, es] = total_emission_i[l-begin_l, h, s, vf, ft, es] + (ef1 * float(street_length[l]) /1000 / speed_data[l,h] + ef2) * flow_data[l,h] \
                                                                           * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                                          fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es]
                                for ef in range(len(emission_factor)):
                                    ### for exhaust
                                    if emission_factor[ef][0] == 'Passenger Cars' and emission_factor[ef][1] == fuel_type_flag_taxi and \
                                            emission_factor[ef][2] == vehicle_fleet_flag_segment[vf] and emission_factor[ef][3] == emission_standard_flag[es] and \
                                            emission_factor[ef][5] == species_flag[s]:
                                        try:
                                            c_speed = ef_equation_china.speed_correction(speed_data[l,h], species_flag[s], emission_factor[ef][1], emission_standard_flag[es])
                                        except UnboundLocalError:
                                            c_speed = 1
                                        try:
                                            c_ta = ef_equation_china.ta_correction(cal_day.month, species_flag[s], emission_factor[ef][1], vehicle_fleet_flag[vf])
                                        except UnboundLocalError:
                                            c_ta = 1
                                        try:
                                            c_rh = ef_equation_china.rh_correction(cal_day.month, species_flag[s], emission_factor[ef][1])
                                        except UnboundLocalError:
                                            c_rh = 1
                                        if high_altitude:
                                            try:
                                                c_height = ef_equation_china.height_correction(species_flag[s],emission_factor[ef][1],vehicle_fleet_flag[vf])
                                            except UnboundLocalError:
                                                c_height = 1
                                        else:
                                            c_height = 1
                                        # c_speed = c_ta = c_rh = c_height = 1
                                        hef_i[l-begin_l,h,s,vf,ft,es] = ef_equation_china.ef_equation_china(float(emission_factor[ef][6]), c_speed, c_ta, c_rh, c_height)
                                        #print(vehicle_fleet_flag[vf], speed_data[l,h], species_flag[s], fuel_type_flag[ft],emission_standard_flag[es], c_speed)
                                        if species_flag[s] == 'NOx':
                                           hef_back=hef_i[l-begin_l,h,s,vf,ft,es]
                                           hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es]*1.5716-0.0882
                                           if hef_i[l-begin_l,h,s,vf,ft,es] <= 0:
                                              hef_i[l-begin_l,h,s,vf,ft,es] = hef_back
                                        ## reference S. Zhang, Y. Wu, H.liu et al, AE, 2013
                                        if fuel_type_flag[ft] == 'LPG':
                                            if species_flag[s] == 'CO':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 0.46
                                            elif species_flag[s] == 'HC':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 1.19
                                            #elif species_flag[s] == 'NOx':
                                            #    hef[l,h,s,vf,ft,es] = hef[l,h,s,vf,ft,es] * 0.75
                                            elif species_flag[s] == 'PM25' or species_flag[s] == 'PM10':
                                                hef_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * 1
                                            else:
                                                hef_i[l-begin_l,h,s,vf,ft,es] = 0
                                        total_emission_i[l-begin_l,h,s,vf,ft,es] = hef_i[l-begin_l,h,s,vf,ft,es] * \
                                                              flow_data[l,h] * vehicle_fleet_percentage[vehicle_fleet_flag[vf]] * \
                                                              fuel_type_percentage[vf, ft] * emission_standard_percentage[vf, es] * \
                                                              float(street_length[l]) /1000
                                        break
        #    count=[]
        #    displ=[]
        #    displ1=0
        #    for count_i in range(0,size):
        #        if count_i==(size-1):
        #            count.append((len(base_coord)-int(math.floor(len(base_coord)//size))*(size-1))*24*len(species_flag)*len(vehicle_fleet_flag)*len(fuel_type_flag)*len(emission_standard_flag))
        #            displ.append(displ1)
        #            displ1+=(len(base_coord)-int(math.floor(len(base_coord)//size))*(size-1))*24*len(species_flag)*len(vehicle_fleet_flag)*len(fuel_type_flag)*len(emission_standard_flag)
        #        else:
        #            count.append(int(math.floor(len(base_coord)//size))*24*len(species_flag)*len(vehicle_fleet_flag)*len(fuel_type_flag)*len(emission_standard_flag))
        #            displ.append(displ1)
        #            displ1+=int(math.floor(len(base_coord)//size))*24*len(species_flag)*len(vehicle_fleet_flag)*len(fuel_type_flag)*len(emission_standard_flag)
        #
        #    comm.Barrier()
        total_emission_ii=np.nansum(np.nansum(np.nansum(total_emission_i, 5), 4), 3)        
        #    if rank==0:
        #        print count 
        #        print displ
        print(total_emission_ii)
        print('now begin to gather')    
        #   comm.Allgatherv([total_emission_i,MPI.DOUBLE],[total_emission,count,displ,MPI.DOUBLE])
        #comm.Gather(total_emission_ii,total_emission_output,root=0)
        #total_emission = comm.allgather(total_emission_i)
        print('Gather is ok')
        #total_emission = np.vstack(total_emission)
        #comm.Barrier()
        #if rank==0:
            #total_emission = result
            #print('gather is ok, now to sum the results')                
            #total_emission_output = np.nansum(np.nansum(np.nansum(total_emission, 5), 4), 3)
        end_time_cal = time.time()
        print('Calculation time:  %.2f s' % (end_time_cal - start_time_cal))
        ### IO part
        start_time_io=time.time()
        for s_out in range(len(species_flag)):
            save_output1 = np.savetxt(output_path + \
                         str(cal_day.year).zfill(4) + '-' + str(cal_day.month).zfill(2) + '-' + str(cal_day.day).zfill(2) +\
                         '_hour_emission_' + species_flag[s_out] + str('%03d'%rank)  + '.txt', \
                                      total_emission_ii[:,:,s_out], fmt='%6f', delimiter=',',
                                      header='0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23', \
                                      newline='\n')
        #for fuel_type_output in range(len(fuel_type_flag)):
        #    scio.savemat(output_path + \
        #                 str(cal_day.year).zfill(4) + '-' + str(cal_day.month).zfill(2) + '-' + str(
        #        cal_day.day).zfill(2) + '_total_emission_' + fuel_type_flag[fuel_type_output] + '.mat', \
        #             {'total_emission_' + fuel_type_flag[fuel_type_output]: total_emission[:,:,:,:,fuel_type_output,:]})
        #for fuel_type_output in range(len(fuel_type_flag)):
        #    scio.savemat(output_path + \
        #                 str(cal_day.year).zfill(4) + '-' + str(cal_day.month).zfill(2) + '-' + str(
        #        cal_day.day).zfill(2) + '_hef_' + fuel_type_flag[fuel_type_output] + '.mat', \
        #             {'hef_' + fuel_type_flag[fuel_type_output]: hef[:,:,:,:,fuel_type_output,:]})
        end_time_io = time.time()
        print('IO time:  %.2f s' % (end_time_io - start_time_io))
        print('Total time:  %.2f s' % (end_time_io - start_time_cal))
