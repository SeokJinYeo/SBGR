# -*- coding: utf-8 -*-
"""
Created on Mon Nov  6 09:28:51 2023

@author: Jin

1. receive the barcode (warping and edge removed)(maybe cell_id assigned) and position files, modify the position file
2. Fine the intersections
3. Calculating the offsets for every intersections
4. Defining duplicated spots


4. Find the cedible offsets (connections) and non-credible offsets
5. Sorting the connections by duplicate spots in the overlapping region
6. Stitching with minimum spanning tree (MST) algorithm
7. Re make the offset network
8. Get the evaluation metrics' scores
9. Removing duplicate spots in credible-major intersections
10. Removing duplicate sopts in uncred0ble-major and minor intersections
11. Save the result

Spot-based global registration (SBGR) - For paper generation
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from skimage.registration import phase_cross_correlation
from shapely import geometry
import networkx as nx
from copy import deepcopy
import warnings
import time
from tqdm import tqdm
from sklearn.neighbors import KDTree
from multiprocessing import Pool
from functools import partial
import time
import multiprocessing
import warnings
from sklearn.neighbors import KNeighborsClassifier
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from scipy.stats import norm

def get_fov_postion_from_zen_fov(zen_fov, positions):
    sel_pos = positions[positions["index"]==zen_fov]
    return [sel_pos["x"].to_numpy()[0], sel_pos["y"].to_numpy()[0]]

def get_fov_bounds1(fov, starting_positions, pix_size):
    min_x, min_y = get_fov_postion_from_zen_fov(fov,starting_positions) 
    max_x = min_x + 2048
    max_y = min_y + 2048
    
    return [min_x*pix_size, min_y*pix_size, max_x*pix_size, max_y*pix_size]

def get_fov_boxes(fovs, starting_positions, pix_size, pix_reduction=0):
    boxes = []
    for f in fovs:
        if pix_reduction!=0:
            min_x, min_y, max_x, max_y = get_fov_bounds1(f,starting_positions,pix_size)
            min_x, min_y = [min_x+pix_reduction*pix_size, min_y+pix_reduction*pix_size]
            max_x, max_y = [max_x-pix_reduction*pix_size, max_y-pix_reduction*pix_size]
            temp_box = geometry.box(min_x, min_y, max_x, max_y)
        else:    
            temp_box = geometry.box(*get_fov_bounds1(f,starting_positions,pix_size))
        boxes.append(temp_box)
    return boxes

def overlapping_box_spot(spots_s,fov,positions,overlapping_box):
    #working_spots=deepcopy(spots)
    #working_spots=working_spots[working_spots['fov']==fov]
    spots_s=spots_s[spots_s['fov']==fov]
    spots_s=spots_s[(spots_s.global_x>=overlapping_box[0])&(spots_s.global_y>=overlapping_box[1])&(spots_s.global_x<=overlapping_box[2])&(spots_s.global_y<=overlapping_box[3])]
    return spots_s

def gene_list(overlapped_spots):
    genes_count=pd.DataFrame(overlapped_spots['barcode_id'].value_counts())
    genes_count.rename(columns={'barcode_id':'count'},inplace=True)
    return genes_count

def find_intensections2(fovs, positions, pix_size):
    fov_boxes = get_fov_boxes(fovs, positions, pix_size)
    intercesting_fovs = []
    minor_intersecting=[]
    for f in fovs:
        for i, x in enumerate(fov_boxes):
            if f != i:
                if fov_boxes[f].intersects(x):
                    overlap_box=fov_boxes[f].intersection(x)
                    x,y=overlap_box.exterior.coords.xy
                    max_len=max(np.linalg.norm(np.array([x[0],y[0]])-np.array([x[1],y[1]])),np.linalg.norm(np.array([x[1],y[1]])-np.array([x[2],y[2]])))
                    if max_len>100:
                        intersec_fov=sorted([f,i])
                        intercesting_fovs.append(intersec_fov)
                    else:
                        intersec_fov=sorted([f,i])
                        minor_intersecting.append(intersec_fov)
    intersecting_fov_pairs=[]
    minor_pairs=[]
    for i in intercesting_fovs:
        if i not in intersecting_fov_pairs:
            intersecting_fov_pairs.append(i)
    for i in minor_intersecting:
        if i not in minor_pairs:
            minor_pairs.append(i)
    return intersecting_fov_pairs,minor_pairs

def gene_ratio(genes_list,gene_list1,gene_list2):
    ratio=np.zeros(len(genes_list))
    for i in range(len(genes_list)):
        if genes_list.index[i] not in gene_list2.index or genes_list.index[i] not in gene_list1.index:
            ratio[i]=0
        else:
            ratio[i]=gene_list2.loc[genes_list.index[i]][0]/(gene_list2.loc[genes_list.index[i]][0]+gene_list1.loc[genes_list.index[i]][0])
    genes_list.insert(1,'ratio',ratio)
    return genes_list
def filtering(genes_list):
    genes_filtered=genes_list[(genes_list['count']>60) &(genes_list['ratio']>0.4) & (genes_list['ratio']<0.6)]
    return genes_filtered

def position_to_image(overlap_gene,box):
    pix_box=box/0.103
    pix_box=pix_box.astype(int)
    #im=np.zeros((pix_box[2]-pix_box[0]+1,pix_box[3]-pix_box[1]+1))
    im=np.zeros((pix_box[2]-pix_box[0]+10,pix_box[3]-pix_box[1]+10))
    size=overlap_gene.shape[0]
    x=overlap_gene['global_x'].to_numpy()
    y=overlap_gene['global_y'].to_numpy()
    #x=x.astype(int)
    #y=y.astype(int)
    for i in range(size):
        x_pix = int((x[i]-box[0])/0.103)
        y_pix = int((y[i]-box[1])/0.103)
        #x_loc=x[i]-box[0]
        #y_loc=y[i]-box[1]
        im[x_pix,y_pix]=255
    return im 

def minimum_spanning_tree(offset_sorted,spots,positions):
    min_tree=pd.DataFrame({'node' : [], 'parent' : [],'rank' : []})
    overlapping_track=[]
    for i in range(len(positions)):
        min_tree.loc[i]=[i,i,1]
    i=0;e=0
    while e<len(positions)-1:
        fov1=offset_sorted.iloc[i]['fov_pair'][0]
        fov2=offset_sorted.iloc[i]['fov_pair'][1]
        pair_offset=np.array([offset_sorted.iloc[i]['x_offset'],offset_sorted.iloc[i]['y_offset']])
        i=i+1
        x=int(min_tree[(min_tree['node']==fov1)]['parent'])
        y=int(min_tree[(min_tree['node']==fov2)]['parent'])
        x_rank=int(min_tree[(min_tree['node']==x)]['rank'])
        y_rank=int(min_tree[(min_tree['node']==y)]['rank'])
        if x!=y:
            e+=1
            overlapping_track.append([fov1, fov2])
            if x_rank < y_rank:
                shift_ind=np.array(min_tree.loc[min_tree['parent']==x,'node'])
                for shift in shift_ind:
                    spots.loc[spots['fov']==shift,'global_x']-=pair_offset[0]*0.103
                    spots.loc[spots['fov']==shift,'global_y']-=pair_offset[1]*0.103
                    positions.loc[positions['index']==shift,'x']-=pair_offset[0]
                    positions.loc[positions['index']==shift,'y']-=pair_offset[1]
                for j in range(i,len(offset_sorted)):
                    f1=offset_sorted.iloc[j]['fov_pair'][0]
                    f2=offset_sorted.iloc[j]['fov_pair'][1]
                    if f1 in shift_ind:
                        offset_sorted.iloc[j,1]-=pair_offset[0]
                        offset_sorted.iloc[j,2]-=pair_offset[1]
                    elif f2 in shift_ind:
                        offset_sorted.iloc[j,1]+=pair_offset[0]
                        offset_sorted.iloc[j,2]+=pair_offset[1]
                min_tree.loc[min_tree['parent']==x,'parent']=y
                min_tree.loc[min_tree['node']==y,'rank']+=x_rank
            elif y_rank <= x_rank:
                shift_ind=np.array(min_tree.loc[min_tree['parent']==y,'node'])
                for shift in shift_ind:
                    spots.loc[spots['fov']==shift,'global_x']+=pair_offset[0]*0.103
                    spots.loc[spots['fov']==shift,'global_y']+=pair_offset[1]*0.103
                    positions.loc[positions['index']==shift,'x']+=pair_offset[0]
                    positions.loc[positions['index']==shift,'y']+=pair_offset[1]
                for j in range(i,len(offset_sorted)):
                    f1=offset_sorted.iloc[j]['fov_pair'][0]
                    f2=offset_sorted.iloc[j]['fov_pair'][1]
                    if f1 in shift_ind:
                        offset_sorted.iloc[j,1]+=pair_offset[0]
                        offset_sorted.iloc[j,2]+=pair_offset[1]
                    elif f2 in shift_ind:
                        offset_sorted.iloc[j,1]-=pair_offset[0]
                        offset_sorted.iloc[j,2]-=pair_offset[1]
                min_tree.loc[min_tree['parent']==y,'parent']=x
                min_tree.loc[min_tree['node']==x,'rank']+=y_rank
    return spots,positions, overlapping_track
    

def removing_dup_spots_knn(spots,positions,offset,intersection_type,error,pix_size,edge_pix,z_step):
    align_log=pd.DataFrame({'intersection' : [], 'count' : [], 'pp':[]})
    z_dif_count=0
    #major credile intersections
    for i in tqdm(range(len(offset))):
        if intersection_type =="minor":
            fov1=offset[i][0]
            fov2=offset[i][1]
        else:
            fov1=offset.iloc[i]['gene_pair'][0] 
            fov2=offset.iloc[i]['gene_pair'][1]
            pair_offset=np.array([offset.iloc[i]['x_offset'],offset.iloc[i]['y_offset']])
            if abs(pair_offset[0])<1 and abs(pair_offset[1])<1:
                pair_offset=[0,0]
        if intersection_type == "un_cred" or intersection_type == "minor" or intersection_type == "evaluate":
            pair_offset=[0,0]   
        box=get_overlapping_box(fov1,fov2,positions,pix_size,-30)
        box=np.array(box.bounds)
        spots_1=overlapping_box_spot(spots,fov1,positions,box)
        spots_2=overlapping_box_spot(spots,fov2,positions,box)
        spots_1_2=pd.concat([spots_1,spots_2])
        genes_list=gene_list(spots_1_2)
        #genes_three=abundant_gene(spots_20_21)
        gene_list1=gene_list(spots_1)
        gene_list2=gene_list(spots_2)
        genes_list=gene_ratio(genes_list,gene_list1,gene_list2)
        overlapped_gene=genes_list[genes_list['ratio']!=0]
        if intersection_type=="evaluate":
            count,count2,pp=align_offset2_knn(spots,fov1,fov2,box,overlapped_gene,pair_offset,pix_size,error,z_step)
        else:
            spots,count,count2,pp=align_offset_knn(spots,fov1,fov2,box,overlapped_gene,pair_offset,pix_size,error,z_step)
            z_dif_count+=count2
        align_log.loc[i]=[[fov1,fov2],count,pp]
    if intersection_type == "cred":
        print("filtered out by z_dif: %i" %z_dif_count)
    return spots,align_log

def get_overlapping_box(fov_0, fov_1, positions_df, pix_size, edge_removal_pix):
    fov_box_0 = get_fov_boxes([fov_0], positions_df, pix_size, edge_removal_pix)[0]
    fov_box_1 = get_fov_boxes([fov_1], positions_df, pix_size, edge_removal_pix)[0]
    
    return fov_box_0.intersection(fov_box_1)
def spots_count_overlapping(offset,spots,positions,edge_pix,pix_size):
    spot_numb=[]
    for minor in range(len(offset)):
        #fov1=minor_intersects[minor][0]
        #fov2=minor_intersects[minor][1]
        fov1=offset.iloc[minor]['gene_pair'][0] 
        fov2=offset.iloc[minor]['gene_pair'][1]
        #fov2=minor[1]
        box=get_overlapping_box(fov1,fov2,positions,pix_size,-30)
        box=np.array(box.bounds)
        spots_1=overlapping_box_spot(spots,fov1,positions,box)
        spots_2=overlapping_box_spot(spots,fov2,positions,box)
        spot_numb.append(len(spots_1)+len(spots_2))
    return spot_numb

def fov_direction(fov1,fov2,positions):
    fov1_pos=np.array(get_fov_postion_from_zen_fov(fov1,positions))
    fov2_pos=np.array(get_fov_postion_from_zen_fov(fov2,positions))
    vector=fov2_pos-fov1_pos
    if abs(vector[0])>abs(vector[1]):
        if vector[0]<0:
            direction=[-1,0]
        elif vector[0]>0:
            direction=[1,0]
    else:
        if vector[1]<0:
            direction=[0,-1]
        else:
            direction=[0,1]
    return direction

def fovs_network_draw(offset,positions):
    fov_align_graph = nx.DiGraph()
    pos={}  
    edge_label={}   
    for i in range(len(positions)):
        fov_align_graph.add_node(i)
        pos[i]=[get_fov_postion_from_zen_fov(i,positions)[0],get_fov_postion_from_zen_fov(i,positions)[1]]
    for j in range(offset.shape[0]):
        node1=offset.iloc[j]['gene_pair'][0]
        node2=offset.iloc[j]['gene_pair'][1]
        #fov_align_graph.add_edge(node1,node2,weight=10)
        fov_align_graph.add_edge(node1,node2,weight=(offset.iloc[j]['x_offset']**2+offset.iloc[j]['y_offset']**2)**0.5)
        edge_label[node1,node2]=[round(offset.iloc[j]['x_offset'],3),round(offset.iloc[j]['y_offset'],3)]
    weights = nx.get_edge_attributes(fov_align_graph,'weight').values()
    #sub figure
    fig = plt.figure()
    fig,ax = plt.subplots(figsize=(70,70))
    #plt.rcParams['figure.figsize'] = (70, 70)  
    nx.draw(fov_align_graph, pos,width=list(weights),node_color='lightgreen', with_labels=True)
    nx.draw_networkx_edge_labels(fov_align_graph, pos, edge_labels=edge_label)
    plt.show()
    
def fovs_network_draw_color_edge(offset,positions,scale):
    fov_align_graph = nx.Graph()
    pos={}  
    edge_label={}   
    for i in range(len(positions)):
        fov_align_graph.add_node(i)
        pos[i]=[get_fov_postion_from_zen_fov(i,positions)[0],get_fov_postion_from_zen_fov(i,positions)[1]]
    for j in range(offset.shape[0]):
        node1=offset.iloc[j]['fov_pair'][0]
        node2=offset.iloc[j]['fov_pair'][1]
        #fov_align_graph.add_edge(node1,node2,weight=10)
        fov_align_graph.add_edge(node1,node2,weight=(offset.iloc[j]['x_offset']**2+offset.iloc[j]['y_offset']**2)**0.5)
        edge_label[node1,node2]=[round(offset.iloc[j]['x_offset'],3),round(offset.iloc[j]['y_offset'],3)]
    a_netw_edges = fov_align_graph.edges()
    a_netw_weights = [fov_align_graph[source][dest]['weight'] for source, dest in a_netw_edges]
    
    
    #scale weights in range 0-1 before assigning color 
    maxWeight=float(max(a_netw_weights))
    minWeight=float(min(a_netw_weights))
    a_netw_colors = [plt.cm.Blues((weight-minWeight)/(maxWeight-minWeight)) for weight in a_netw_weights]
    plt.ioff()
    #multiply all tuples in color list by scale factor
    #colors_unscaled=[tuple(map(lambda x: (maxWeight-minWeight)*x+minWeight-2, y)) for y in a_netw_colors]
    # um scale
    if scale == 'um':
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight)*0.103, y)) for y in a_netw_colors]
    elif scale == 'nm':
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight)*103, y)) for y in a_netw_colors]
    else:
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight), y)) for y in a_netw_colors]
    #generate a 'dummy' heatmap using the edgeColors as substrate for colormap
    heatmap = plt.pcolor(colors_unscaled,cmap=plt.cm.Blues)
    
    #re-enable plotting
    plt.ion()
    #sub figure
    fig = plt.figure()
    fig,ax = plt.subplots(figsize=(10,10))
    #fig,ax = plt.subplots()
    #plt.rcParams['figure.figsize'] = (70, 70)  
    #nx.draw(fov_align_graph, pos,width=list(weights),node_color='lightgreen', edge_color=a_netw_colors,with_labels=True, ax=ax)
    nx.draw(fov_align_graph, pos,width=10,node_color='lightgreen', edge_color=a_netw_colors,with_labels=False, ax=ax, node_size = 100)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect("equal")
    cbar = plt.colorbar(heatmap)
    cbar.ax.tick_params(labelsize=15)
    #cbar.ax.set_ylabel('edge weight',labelpad=100,rotation=270, fontsize=100)
    #nx.draw_networkx_edge_labels(fov_align_graph, pos, edge_labels=edge_label)
    plt.show()
    return maxWeight, minWeight, colors_unscaled
    
def fovs_network_draw_color_edge2(offset,positions, maxWeight, minWeight,scale):
    fov_align_graph = nx.Graph()
    pos={}  
    edge_label={}   
    for i in range(len(positions)):
        fov_align_graph.add_node(i)
        pos[i]=[get_fov_postion_from_zen_fov(i,positions)[0],get_fov_postion_from_zen_fov(i,positions)[1]]
    for j in range(offset.shape[0]):
        node1=offset.iloc[j]['fov_pair'][0]
        node2=offset.iloc[j]['fov_pair'][1]
        #fov_align_graph.add_edge(node1,node2,weight=10)
        fov_align_graph.add_edge(node1,node2,weight=(offset.iloc[j]['x_offset']**2+offset.iloc[j]['y_offset']**2)**0.5)
        edge_label[node1,node2]=[round(offset.iloc[j]['x_offset'],3),round(offset.iloc[j]['y_offset'],3)]
    a_netw_edges = fov_align_graph.edges()
    a_netw_weights = [fov_align_graph[source][dest]['weight'] for source, dest in a_netw_edges]
    a_netw_colors = [plt.cm.Blues((weight-minWeight)/(maxWeight-minWeight)) for weight in a_netw_weights]
    if scale == 'um':
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight)*0.103, y)) for y in a_netw_colors]
    elif scale == 'nm':
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight)*103, y)) for y in a_netw_colors]
    else:
        colors_unscaled=[tuple(map(lambda x: ((maxWeight-minWeight)*x+minWeight), y)) for y in a_netw_colors]
    #colors_unscaled=[tuple(map(lambda x: (maxWeight-minWeight)*x+minWeight-2, y)) for y in a_netw_colors]
    plt.ioff()
    #multiply all tuples in color list by scale factor
    #generate a 'dummy' heatmap using the edgeColors as substrate for colormap
    heatmap = plt.pcolor(colors_unscaled,cmap=plt.cm.Blues)
    
    #re-enable plotting
    plt.ion()
    #sub figure
    fig = plt.figure()
    fig,ax = plt.subplots(figsize=(10,10))
    #fig,ax = plt.subplots()
    #plt.rcParams['figure.figsize'] = (70, 70)  
    #nx.draw(fov_align_graph, pos,width=list(weights),node_color='lightgreen', edge_color=a_netw_colors,with_labels=True, ax=ax)
    nx.draw(fov_align_graph, pos,width=10,node_color='lightgreen', edge_color = a_netw_colors,with_labels=False, ax=ax, node_size = 100)
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.set_aspect("equal")
    cbar = plt.colorbar(heatmap)
    cbar.ax.tick_params(labelsize=15)
    #cbar.ax.set_ylabel('edge weight',labelpad=100,rotation=270, fontsize=100)
    #nx.draw_networkx_edge_labels(fov_align_graph, pos, edge_labels=edge_label)
    plt.show()

def fov_coord(loaded_positions_df):
    x = set(list(loaded_positions_df['x']))
    x = sorted(x)
    x_coord = {}
    for i in range(len(x)):
        x_coord[x[i]] = i
        
    y = set(list(loaded_positions_df['y']))
    y = sorted(y)
    y_coord = {}
    for i in range(len(y)):
        y_coord[y[i]] = i
        
    fov_coord = {}
    for i in range(len(loaded_positions_df)):
        tmp_x_coord = x_coord[loaded_positions_df.iloc[i]['x']]
        tmp_y_coord = y_coord[loaded_positions_df.iloc[i]['y']]
        fov_coord[loaded_positions_df.iloc[i]['index']] = np.array([tmp_x_coord,tmp_y_coord])
    return fov_coord

def align_offset_knn(spots,fov1,fov2,overlapping_box,genes_list,offset,pix_size,error,z_step):
    count=0
    count2=0
    pp=[]
    spots2=spots[(spots.fov==fov2)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots1=spots[(spots.fov==fov1)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots2['global_x']+=offset[0]*pix_size
    spots2['global_y']+=offset[1]*pix_size
    
    for i in range(len(genes_list)):
        tmp_spot = spots1[spots1['barcode_id']==genes_list.index[i]]
        tmp_spot2 = spots2[spots2['barcode_id']==genes_list.index[i]]

        tmp_spots_12 = pd.concat([tmp_spot,tmp_spot2])
        index_list = list(tmp_spots_12.index)
        
        x=np.array(tmp_spots_12['global_x'])
        y=np.array(tmp_spots_12['global_y'])
        z = np.zeros((len(tmp_spots_12),2))
        z[:,0]=x
        z[:,1]=y
        
        tree = KDTree(z)
        dist, ind = tree.query(z,k=len(tmp_spots_12))  
        dist_1 = dist[:len(tmp_spot)]
        ind_1 = ind[:len(tmp_spot)]
        ind_1_t = ind_1>=len(tmp_spot)
        
        for j in range(len(tmp_spot)):
            if dist_1[j,:][ind_1_t[j,:]][0] <= error: #dup spot
                spot1_ind = index_list[j]
                spot2_closet_ind = index_list[ind_1[j,:][ind_1_t[j,:]][0]]
                spot_pos1 = np.array([tmp_spot.loc[spot1_ind,'global_x'],tmp_spot.loc[spot1_ind,'global_y']])
                spot_pos2 = np.array([tmp_spot2.loc[spot2_closet_ind,'global_x'],tmp_spot2.loc[spot2_closet_ind,'global_y']])
                spot_pos_final = np.array([(spot_pos2[0]+spot_pos1[0])/2,(spot_pos2[1]+spot_pos1[1])/2])
                z_dif=z_step*abs(tmp_spot.loc[spot1_ind,'global_z']-tmp_spot2.loc[spot2_closet_ind,'global_z'])
                if z_dif<3:
                    spot_pos_final_z = (tmp_spot.loc[spot1_ind,'global_z']+tmp_spot2.loc[spot2_closet_ind,'global_z'])/2
                    spots.loc[spot1_ind,'barcode_id']=-1
                    spots.loc[spot2_closet_ind,'global_x']=spot_pos_final[0]
                    spots.loc[spot2_closet_ind,'global_y']=spot_pos_final[1]
                    spots.loc[spot2_closet_ind,'global_z']=spot_pos_final_z
                    pp.append(spot1_ind)
                    #print(fov1_box_index[0])
                    count+=1
                else:
                    count2+=1
    return spots,count,count2,pp

def align_offset2_knn(spots,fov1,fov2,overlapping_box,genes_list,offset,pix_size,error,z_step):
    count=0
    count2=0
    pp=[]
    spots2=spots[(spots.fov==fov2)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots1=spots[(spots.fov==fov1)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
	

    spots2['global_x']+=offset[0]*pix_size
    spots2['global_y']+=offset[1]*pix_size
    
    for i in range(len(genes_list)):
        tmp_spot = spots1[spots1['barcode_id']==genes_list.index[i]]
        tmp_spot2 = spots2[spots2['barcode_id']==genes_list.index[i]]

        tmp_spots_12 = pd.concat([tmp_spot,tmp_spot2])
        index_list = list(tmp_spots_12.index)
        
        x=np.array(tmp_spots_12['global_x'])
        y=np.array(tmp_spots_12['global_y'])
        z = np.zeros((len(tmp_spots_12),2))
        z[:,0]=x
        z[:,1]=y
        
        tree = KDTree(z,leaf_size=40)
        dist, ind = tree.query(z,k=len(tmp_spots_12),dualtree=True)  
        dist_1 = dist[:len(tmp_spot)]
        ind_1 = ind[:len(tmp_spot)]
        ind_1_t = ind_1>=len(tmp_spot)
        
        for j in range(len(tmp_spot)):
            if dist_1[j,:][ind_1_t[j,:]][0] <= error: #dup spot
                spot1_ind = index_list[j]
                spot2_closet_ind = index_list[ind_1[j,:][ind_1_t[j,:]][0]]
                z_dif=z_step*abs(tmp_spot.loc[spot1_ind,'global_z']-tmp_spot2.loc[spot2_closet_ind,'global_z'])
                if z_dif<3:
                    count+=1
                    pp.append(spot1_ind)
                else:
                    count2+=1
    return count,count2,pp

    

# Calculating Offset

def offset_calc(working_spots,positions, pix_size, fovs):
    i=0
    fov1=fovs[0]
    fov2=fovs[1]
    offset_df=pd.DataFrame({'fov_pair' : [], 'x_offset' : [],'y_offset' : [],'tot_spot':[]})
    box=get_overlapping_box(fov1,fov2,positions,pix_size,-30)
    #box=get_overlapping_box(fov1,fov2,positions,pix_size,0)
    box=np.array(box.bounds)
    spots_1=overlapping_box_spot(working_spots,fov1,positions,box)
    spots_2=overlapping_box_spot(working_spots,fov2,positions,box)
    im1_1=position_to_image(spots_1,box)
    im2_1=position_to_image(spots_2,box)
    warnings.filterwarnings(action='ignore')
    offset1=phase_cross_correlation(im1_1,im2_1,upsample_factor=100)
    offset_df.loc[i]=[fovs,offset1[0][0],offset1[0][1],len(spots_1)+len(spots_2)]    
    return offset_df

def offset_calc_parallel(spots,positions,intersecting_fov_pairs, pix_size):
    offset_df=pd.DataFrame({'fov_pair' : [], 'x_offset' : [],'y_offset' : [],'tot_spot':[]})
    working_spots = spots[['fov','global_x','global_y','global_z','barcode_id','id']]
    
    pool = Pool(processes=multiprocessing.cpu_count()) 
    custfunct = partial(offset_calc,working_spots,positions,pix_size)
    fov_pair_list=[]
    for i in range(len(intersecting_fov_pairs)):
        fov_pair_list.append([intersecting_fov_pairs[i][0],intersecting_fov_pairs[i][1]])

    result = pool.map_async(custfunct, fov_pair_list, chunksize=100)
    for value in result.get():
        offset_df=pd.concat([offset_df,value])
    pool.close()
    return offset_df

def spots_to_pos(spots):
    ids = spots['id']
    x = spots['global_x']
    y = spots['global_y']
    df = pd.DataFrame({'id' : ids,'x' : x, 'y' : y})
    return df
    
    

def spot_dist2_knn2(spots,fov1,fov2,overlapping_box,genes_list,offset,pix_size):
    spots2=spots[(spots.fov==fov2)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots1=spots[(spots.fov==fov1)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots2['global_x']+=offset[0]*pix_size
    spots2['global_y']+=offset[1]*pix_size
    dup_df=pd.DataFrame({'id1' : [], 'id2' : [],'dist' : [], 'fov_pair':[]})
    for i in range(len(genes_list)):
        
        tmp_spot = spots1[spots1['barcode_id']==genes_list.index[i]]
        tmp_spot2 = spots2[spots2['barcode_id']==genes_list.index[i]]
        
        pos1 = spots_to_pos(tmp_spot)
        pos2 = spots_to_pos(tmp_spot2)

        if len(pos2)<len(pos1):
            [nn_distance, nn_spot_ids] = knn_imputation(pos1,pos2,1)
            tmp_id = np.array(pos1['id'])
            nn_df_id =[int(tmp_id[i]) for i in nn_spot_ids]
            if i==0:
                dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos2['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id)})
            else:
                tmp_dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos2['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id)})
                dup_df = pd.concat([dup_df,tmp_dup_df])
            
        else:
            [nn_distance, nn_spot_ids] = knn_imputation(pos2,pos1,1)
            tmp_id = np.array(pos2['id'])
            nn_df_id =[int(tmp_id[i]) for i in nn_spot_ids]
            if i==0:
                dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos1['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id)})
            else:
                tmp_dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos1['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id)})
                dup_df = pd.concat([dup_df,tmp_dup_df])
    return dup_df

def spot_dist2_knn3(spots,fov1,fov2,overlapping_box,genes_list,offset,pix_size):
    spots2=spots[(spots.fov==fov2)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots1=spots[(spots.fov==fov1)&(spots.global_x>=overlapping_box[0])&(spots.global_y>=overlapping_box[1])&(spots.global_x<=overlapping_box[2])&(spots.global_y<=overlapping_box[3])]
    spots2['global_x']+=offset[0]*pix_size
    spots2['global_y']+=offset[1]*pix_size
    dup_df=pd.DataFrame({'id1' : [], 'id2' : [],'dist' : [], 'fov_pair':[],'x1':[],'x2':[],'y1':[],'y2':[]})
    for i in range(len(genes_list)):
        
        tmp_spot = spots1[spots1['barcode_id']==genes_list.index[i]]
        tmp_spot2 = spots2[spots2['barcode_id']==genes_list.index[i]]
        
        pos1 = spots_to_pos(tmp_spot)
        pos2 = spots_to_pos(tmp_spot2)
        if len(pos2)<len(pos1):
            [nn_distance, nn_spot_ids] = knn_imputation(pos1,pos2,1)
            tmp_id = np.array(pos1['id'])
            x1 = np.array(pos1['x'])
            y1 = np.array(pos1['y'])
            nn_df_id =[int(tmp_id[i]) for i in nn_spot_ids]
            nn_x = [int(x1[i]) for i in nn_spot_ids]
            nn_y = [int(y1[i]) for i in nn_spot_ids]
            if i==0:
                dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos2['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id),'x1':nn_x,'x2':np.array(pos2['x']),'y1':nn_y,'y2':np.array(pos2['y'])})
            else:
                tmp_dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos2['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id),'x1':nn_x,'x2':np.array(pos2['x']),'y1':nn_y,'y2':np.array(pos2['y'])})
                dup_df = pd.concat([dup_df,tmp_dup_df])
            
        else:
            [nn_distance, nn_spot_ids] = knn_imputation(pos2,pos1,1)
            tmp_id = np.array(pos2['id'])
            nn_df_id =[int(tmp_id[i]) for i in nn_spot_ids]
            x1 = np.array(pos2['x'])
            y1 = np.array(pos2['y'])
            nn_x = [int(x1[i]) for i in nn_spot_ids]
            nn_y = [int(y1[i]) for i in nn_spot_ids]
            if i==0:
                dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos1['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id),'x1':nn_x,'x2':np.array(pos1['x']),'y1':nn_y,'y2':np.array(pos1['y'])})
            else:
                tmp_dup_df = pd.DataFrame({'id1' : nn_df_id,'id2' : np.array(pos1['id']), 'dist' : list(nn_distance[:,0]),'fov_pair':[[fov1,fov2]]*len(nn_df_id),'x1':nn_x,'x2':np.array(pos1['x']),'y1':nn_y,'y2':np.array(pos1['y'])})
                dup_df = pd.concat([dup_df,tmp_dup_df])
    return dup_df


def two_fov_dist(spots,positions,pix_size,vals):
    dup_df=pd.DataFrame({'id1' : [], 'id2' : [],'dist' : [], 'fov_pair':[],'x1':[],'x2':[],'y1':[],'y2':[]})
    fov1=vals[0]
    fov2=vals[1]
    pair_offset=vals[2]
    box=get_overlapping_box(fov1,fov2,positions,pix_size,0)
    box=np.array(box.bounds)
    spots_1=overlapping_box_spot(spots,fov1,positions,box)
    spots_2=overlapping_box_spot(spots,fov2,positions,box)
    spots_1_2=pd.concat([spots_1,spots_2])
    genes_list=gene_list(spots_1_2)
    #genes_three=abundant_gene(spots_20_21)
    gene_list1=gene_list(spots_1)
    gene_list2=gene_list(spots_2)
    genes_list=gene_ratio(genes_list,gene_list1,gene_list2)
    overlapped_gene=genes_list[genes_list['ratio']!=0]
    dup_df = spot_dist2_knn3(spots,fov1,fov2,box,overlapped_gene,pair_offset,pix_size)
    return dup_df
    

def distance_dist2_parallel(spots,positions,offset,pix_size):
    warnings.filterwarnings(action='ignore')
    working_spots = spots[['fov','global_x','global_y','global_z','barcode_id','id']]
    dup_dists_df=pd.DataFrame({'id1' : [], 'id2' : [],'dist' : [], 'fov_pair':[],'x1':[],'x2':[],'y1':[],'y2':[]})
    #major credile intersections
    for i in tqdm(range(len(offset))):
        fov1=offset.iloc[i]['fov_pair'][0] 
        fov2=offset.iloc[i]['fov_pair'][1]
        pair_offset=np.array([offset.iloc[i]['x_offset'],offset.iloc[i]['y_offset']])
        tmp_df = two_fov_dist(working_spots,positions,pix_size,[fov1,fov2,pair_offset])
        dup_dists_df = pd.concat([dup_dists_df,tmp_df])
    return dup_dists_df


def knn_imputation (reference_data, target_data, k_nn, params={}):
    neigh = KNeighborsClassifier(n_neighbors=k_nn, n_jobs=-1,
                                 metric='euclidean', algorithm='brute')
    neigh.fit(reference_data[["x","y"]], np.array(reference_data["id"]))
    nn_distance, nn_cell_ids = neigh.kneighbors(target_data[["x","y"]])
    return [nn_distance, nn_cell_ids]

def spot_dist_calc(spots,id1,id2):
    df1 = spots.loc[[id1]]
    df2 = spots.loc[[id2]]
    return ((float(df1['global_x'])-float(df2['global_x']))**2 + (float(df1['global_y'])-float(df2['global_y']))**2)**0.5

def spot_cnt(working_spots,positions, pix_size,tmp_df):
    fovs = tmp_df['fov_pair'][0]
    offset = np.array([tmp_df['x_offset'][0],tmp_df['y_offset'][0]])
    dup_spots = int(tmp_df['dup_spots'])
    fov1=fovs[0]
    fov2=fovs[1]
    sel_pos = positions[positions["index"]==fov1]
    min_x, min_y =  [sel_pos["x"].to_numpy()[0], sel_pos["y"].to_numpy()[0]]
    max_x = min_x + 2048
    max_y = min_y + 2048
    box1 = geometry.box(*[(min_x+30)*pix_size, (min_y+30)*pix_size, (max_x-30)*pix_size, (max_y-30)*pix_size])
    
    sel_pos = positions[positions["index"]==fov2]
    min_x, min_y =  [sel_pos["x"].to_numpy()[0]+offset[0], sel_pos["y"].to_numpy()[0]+offset[1]]
    max_x = min_x + 2048
    max_y = min_y + 2048
    box2 = geometry.box(*[(min_x+30)*pix_size, (min_y+30)*pix_size, (max_x-30)*pix_size, (max_y-30)*pix_size])
    
    box=np.array(box1.intersection(box2).bounds)
    if len(box)!=0:
        spots_1=overlapping_box_spot(working_spots,fov1,positions,box)
        
        spots2=working_spots[working_spots['fov']==fov2]
        spots2['global_x']+=offset[0]*pix_size
        spots2['global_y']+=offset[1]*pix_size
        spots_2=overlapping_box_spot(spots2,fov2,positions,box)
        spot_numb = len(spots_1)+len(spots_2)
        
        #ncc
        im1_2=position_to_image(spots_1,box)
        im2_2=position_to_image(spots_2,box)
        im1b=im1_2-np.mean(im1_2)
        im2b=im2_2-np.mean(im2_2)
        denominator=(np.sum(im1b*im1b)**0.5)*(np.sum(im2b*im2b)**0.5)
        if denominator==0:
            norm_cross=float("NaN")
        else:
            norm_cross=np.sum(im1b*im2b)/denominator
            
        if spot_numb!=0 and dup_spots!='nan':
            dup_ratio = 2*dup_spots/spot_numb
        else:
            dup_ratio = 0
            
        tmp_df['ncc'] = norm_cross
        tmp_df['dup_ratio'] = dup_ratio
        tmp_df['tot_spots_after_shift'] = spot_numb
    else:
        tmp_df['ncc'] = 0
        tmp_df['dup_ratio'] = 0
        tmp_df['tot_spots_after_shift'] = 0
    return tmp_df
    
def spot_cnt_parallel(spots,positions, pix_size ,offset_df):
    working_spots = deepcopy(spots)
    pool = Pool(processes=multiprocessing.cpu_count()) 
    custfunct = partial(spot_cnt,working_spots,positions,pix_size)
    job_list=[]
    for i in range(len(offset_df)):
        tmp_df = offset_df.iloc[[i]]
        job_list.append(tmp_df)

    result = pool.map_async(custfunct, job_list, chunksize=100)
    new_df = offset_df.iloc[[i]]
    for value in result.get():
        new_df=pd.concat([new_df,value])
    pool.close()
    new_df = new_df.iloc[1: , :]
        
    return new_df

def spot_cnt_non_parallel(spots,positions, pix_size ,offset_dd):
    working_spots = deepcopy(spots)
    job_list=[]
    new_df = offset_dd.iloc[[0]]
    for i in range(len(offset_dd)):
        print(i)
        tmp_df = offset_dd.iloc[[i]]
        tmp_df = spot_cnt(working_spots,positions,pix_size,tmp_df)
        new_df=pd.concat([new_df,tmp_df])

    new_df = new_df.iloc[1: , :]  
    return new_df

def overlapping_eval(spots,positions,fovs,pix_size):
    #Calculating offsets
    start = time.time()
    offset_df=offset_calc_parallel(spots,positions,fovs,pix_size)
    end = time.time()
    print(f"{end - start:.5f} sec for calculating offsets")
    
    # Aggregating spot pair-wise distances
    start = time.time()
    dup_distances_df = distance_dist2_parallel(spots,positions,offset_df,pix_size)
    end = time.time()
    print(f"{end - start:.5f} sec for aggregating spot pair distances")
    
    dup_distances_df.index = range(len(dup_distances_df))
    min_indices = dup_distances_df.groupby('id1')['dist'].idxmin()
    final_dup_dist_df = dup_distances_df.loc[min_indices]         
            
    
    
    dist_array = np.array(final_dup_dist_df['dist'])
    dist_log = np.log10(dist_array)

    # Create sample data

    # Fit the GMM model to the data
    gmm = GaussianMixture(n_components=2)
    gmm.fit(dist_log.reshape(-1, 1))
    mean = gmm.means_  
    print("gaussian fit menas",mean[0][0],mean[1][0])
    covs  = gmm.covariances_
    weights = gmm.weights_
    x_axis = np.arange(-3,3, 0.01)
    y_axis0 = norm.pdf(x_axis, float(mean[0][0]), np.sqrt(float(covs[0][0][0])))*weights[0] # 1st gaussian
    y_axis1 = norm.pdf(x_axis, float(mean[1][0]), np.sqrt(float(covs[1][0][0])))*weights[1] # 2nd gaussian
    plt.rcParams.update({'figure.dpi':500})
    fig,ax = plt.subplots()
    plt.xlim(-2.5, 2.5)
    plt.hist(dist_log, color='black',density = True, bins=100)
    plt.plot(x_axis, y_axis0, lw=3, c='C0')
    plt.plot(x_axis, y_axis1, lw=3, c='C1')
    plt.xlabel("log(distance)", fontsize=20)
    plt.ylabel("Density", fontsize=20)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.xlim(-2.5, 2.5)
    plt.hist(dist_log, color='black', bins=100)
    plt.xlabel("log(distance)", fontsize=20)
    plt.ylabel("Counts", fontsize=20)
    plt.show()
    
    # Get the cluster assignments for each data point
    clusters = gmm.predict(dist_log.reshape(-1, 1))
    
    dist_df = pd.DataFrame({'dist':dist_log, 'label' : clusters})
    dist_0 = dist_df[dist_df['label']==0]
    dist_1 = dist_df[dist_df['label']==1]
    final_dup_dist_df['clusters']=clusters
    
    dist_df2 = pd.DataFrame({'dist':dist_array, 'label' : clusters})
    dist_0 = dist_df2[dist_df2['label']==0]
    dist_1 = dist_df2[dist_df2['label']==1]
    if np.mean(dist_0['dist']) > np.mean(dist_1['dist']):
        error1 = np.max(dist_1['dist'])
        error1_eval = np.mean(dist_1['dist'])+2*np.std(dist_1['dist'])
        identified_dup_spots = final_dup_dist_df[final_dup_dist_df['clusters']==1]
        if error1_eval > 1:
            print(error1_eval)
            print("clustering error")
    else:
        error1 = np.max(dist_0['dist'])
        error1_eval = np.mean(dist_0['dist'])+2*np.std(dist_0['dist'])
        identified_dup_spots = final_dup_dist_df[final_dup_dist_df['clusters']==0]
        if error1_eval > 1:
            print(error1_eval)
            print("clustering error") 
    print("threshold : ", error1)
    print("means : ", np.mean(dist_0['dist']), np.mean(dist_1['dist']))
    print("stds : ", np.std(dist_0['dist']), np.std(dist_1['dist']))
    # Calculating 
    ori_dist = []
    for i in tqdm(range(len(identified_dup_spots))):
        tmp_df = identified_dup_spots.iloc[[i]]
        ori_dist.append(spot_dist_calc(spots,int(tmp_df['id1']),int(tmp_df['id2'])))
        
    identified_dup_spots['ori_dist'] = ori_dist
    identified_dup_spots['fov_pair_str'] = [str(i) for i in identified_dup_spots['fov_pair']]
    
    d_err = []
    dup_spots = []
    for tmp_pair in a:
        identified_dup_spots_ss = identified_dup_spots[identified_dup_spots['fov_pair_str']==str(tmp_pair)]
        d_err.append(np.average(identified_dup_spots_ss['ori_dist']))
        dup_spots.append(len(identified_dup_spots_ss))
        
    offset_df['d_err'] = d_err
    offset_df['dup_spots'] = dup_spots
    
    offset_df = spot_cnt_parallel(spots, positions, pix_size,offset_df)
    
    
    translations = (np.array(offset_df['x_offset']*pix_size)**2 + np.array(offset_df['y_offset']*pix_size)**2)**0.5
    offset_df['translations'] = translations
    offset_non_cred_enlarge = offset_df[offset_df['dup_spots']<11]
    
    fig,ax = plt.subplots()
    plt.hist(offset_df['dup_ratio'], color='black', bins=10)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.hist(offset_df['translations'], color='black', bins=20)
    plt.ylabel("Coutns", fontsize=20)
    plt.xlabel("Translations[um]", fontsize=20)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.hist(offset_df['translations'], color='black', bins=50)
    plt.ylabel("Coutns", fontsize=20)
    plt.xlabel("Translations[um]", fontsize=20)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.hist(offset_df['translations'], color='black', bins=100)
    plt.ylabel("Coutns", fontsize=20)
    plt.xlabel("Translations[um]", fontsize=20)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_df['dup_spots'],offset_df['dup_ratio'], s= 5,c = offset_df['translations'],cmap='gist_rainbow',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("Ratio", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(np.arange(0,1,0.2),fontsize=18)
    
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_df['dup_spots'],offset_df['dup_ratio'], s= 5,c = offset_df['translations'],cmap='jet',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("Ratio", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(np.arange(0,1,0.2),fontsize=18)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_non_cred_enlarge['dup_spots'],offset_non_cred_enlarge['dup_ratio'], s= 15,c = offset_non_cred_enlarge['translations'],cmap='gist_rainbow',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("Ratio", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(np.arange(0,1,0.2),fontsize=18)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_non_cred_enlarge['dup_spots'],offset_non_cred_enlarge['dup_ratio'], s= 5,c = offset_non_cred_enlarge['translations'],cmap='jet',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("Ratio", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(fontsize=18)
    plt.yticks(np.arange(0,1,0.2),fontsize=18)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_df['dup_spots'],offset_df['ncc'], s= 5,c = offset_df['translations'],cmap='gist_rainbow',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("NCC", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_df['dup_spots'],offset_df['ncc'], s= 5,c = offset_df['translations'],cmap='jet',vmin=0,vmax=100)
    plt.colorbar()
    plt.ylabel("NCC", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(fontsize=18)
    plt.show()
    
    cred_or_not = []
    for i in range(len(offset_df)):
        tmp_df = offset_df.iloc[[i]]
        if float(tmp_df['dup_ratio'])>0.2 and int(tmp_df['dup_spots'])>2:
            cred_or_not.append(1)
        else:
            cred_or_not.append(0)
            
    offset_df['cred'] = cred_or_not
    
    offset_cred = offset_df[offset_df['cred']==1]
    offset_non_cred = offset_df[offset_df['cred']==0]
    
    
    
    fig,ax = plt.subplots()
    plt.scatter(offset_cred['dup_spots'],offset_cred['dup_ratio'], s= 5,c = 'blue',label = 'Credible')
    plt.scatter(offset_non_cred['dup_spots'],offset_non_cred['dup_ratio'], s= 5,c = 'black', label = 'Non-credible')
    plt.ylabel("Ratio", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(np.arange(0,1,0.2),fontsize=18)
    plt.colorbar()
    #plt.legend(fontsize=10)
    plt.show()
    
    fig,ax = plt.subplots()
    plt.scatter(offset_cred['dup_spots'],offset_cred['ncc'], s= 5,c = 'blue',label = 'Credible')
    plt.scatter(offset_non_cred['dup_spots'],offset_non_cred['ncc'], s= 5,c = 'black', label = 'Non-credible')
    plt.ylabel("NCC", fontsize=20)
    plt.xlabel("Duplicated spot pairs", fontsize=20)
    plt.xticks(np.arange(0,1300,500),fontsize=18)
    plt.yticks(fontsize=18)
    plt.colorbar()
    #plt.legend(fontsize=10)
    plt.show()
    
    return offset_df,offset_cred,offset_non_cred,final_dup_dist_df,identified_dup_spots, error1
    


if __name__ == '__main__':
#parameters set-up
    save_folder = r'D:\Bee-brain\20220630\cs3'
    pix_size=0.103 #um
    edge_pix=30 #removed edge pixel size
    warped=0 #unwarped spots
    z_step = 1.5 #um
    
    #receive potision.csv file and trim
    positions = pd.read_csv('//core-server.igb.illinois.edu/groups/hshan/sy44/fine_alignment_iteration_datasets/20220630_cs3/positions.csv', header=None, names=["x", "y"])
    #positions = pd.read_csv('D:/U2OS/cs2_strip_mfish/merlin_output/spots/positions.csv', header=None, names=["x", "y"])
    positions["x"] = positions["x"]/pix_size 
    positions["y"] = positions["y"]/pix_size
    fov_info=np.zeros((len(positions),2))
    positions.insert(0, "index", np.arange(0,len(positions)), True)
    
    #receiving the spots file
    spots=pd.read_csv('//core-server.igb.illinois.edu/groups/hshan/sy44/fine_alignment_iteration_datasets/20220630_cs3/barcodes_20230501_warped.csv')
    spots['id'] = spots.index
    fovs=range(len(positions))
    
    #find the major(a) and minor(b) intersections
    a,b=find_intensections2(fovs,positions,pix_size)
    
    offset_df,offset_cred,offset_non_cred,final_dup_dist_df,identified_dup_spots,error1 = overlapping_eval(spots,positions,a,pix_size)
    offset_df_ori  = deepcopy(offset_df)
    offset_df.to_csv(save_folder + '/offsets_0rnd.csv')
    identified_dup_spots.to_csv(save_folder + '/dup_spots_0rnd.csv')
    
    #network generation and estimation
    offset_df_before_est = deepcopy(offset_df)
    offset_df.index = range(len(offset_df))
    warped = 0 # Is it already stitched by other stitching method
    if warped==1:
        offset_df.loc[(offset_df['cred']==0),'x_offset']=0
        offset_df.loc[(offset_df['cred']==0),'y_offset']=0
    else:
        fov_cord = fov_coord(positions)
        
        horizontal_offsets=[]
        vertical_offsets = {}
        for i in range(len(offset_df)):
            pair=offset_df.iloc[i]['fov_pair']
            if offset_df.iloc[i]['cred']==1:
                direction=fov_cord[pair[1]]-fov_cord[pair[0]]
                if str(direction) == str(np.array([0,1])):
                    if fov_cord[pair[0]][1] not in vertical_offsets:
                        vertical_offsets[fov_cord[pair[0]][1]] = [np.array([offset_df.iloc[i]['x_offset'],offset_df.iloc[i]['y_offset']])]
                    else:
                        vertical_offsets[fov_cord[pair[0]][1]].append(np.array([offset_df.iloc[i]['x_offset'],offset_df.iloc[i]['y_offset']]))
                else:
                    if str(direction) == str(np.array([1,0])):
                        horizontal_offsets.append(np.array([offset_df.iloc[i]['x_offset'],offset_df.iloc[i]['y_offset']]))
                    if str(direction) == str(np.array([-1,0])):
                        horizontal_offsets.append(np.array([-offset_df.iloc[i]['x_offset'],-offset_df.iloc[i]['y_offset']]))
        
                    
        mean_vertical_offset = {}
        for key in vertical_offsets:
            mean_vertical_offset[key]=np.mean(vertical_offsets[key],axis=0)
        mean_horizontal_offset = np.mean(horizontal_offsets,axis=0)
    
        for i,row in offset_df.iterrows():
            pair=offset_df.loc[i,'fov_pair']
            if offset_df.loc[i,'cred']==0:
                direction=fov_cord[pair[1]]-fov_cord[pair[0]]
                if str(direction) == str(np.array([0,1])):
                    offset_df.loc[i,'x_offset']=mean_vertical_offset[fov_cord[pair[0]][1]][0]
                    offset_df.loc[i,'y_offset']=mean_vertical_offset[fov_cord[pair[0]][1]][1]
                else:
                    if str(direction) == str(np.array([1,0])):
                        offset_df.loc[i,'x_offset']=mean_horizontal_offset[0]
                        offset_df.loc[i,'y_offset']=mean_horizontal_offset[1]
                    if str(direction) == str(np.array([-1,0])):
                        offset_df.loc[i,'x_offset']=-mean_horizontal_offset[0]
                        offset_df.loc[i,'y_offset']=-mean_horizontal_offset[1]
                        
    maxweight,minweight,scale = fovs_network_draw_color_edge(offset_cred,positions,'um')
    fovs_network_draw_color_edge2(offset_df,positions, maxweight,minweight,'um')
    
    offset_sorted=offset_df.sort_values('dup_spots',ascending=False) # weighted by duplicate spots
    #offset_sorted=offset.sort_values('duplicate',ascending=True) # reverse weighted by duplicate spots
    spots,positions,overlapping_track = minimum_spanning_tree(offset_sorted,spots,positions)
    
    offset_df_2,offset_cred_2,offset_non_cred_2,final_dup_dist_df_2,identified_dup_spots_2,error2 = overlapping_eval(spots,positions,a,pix_size)
    offset_df_2['cred_bef'] = list(offset_df_ori['cred'])
    
    offset_df_2_cred = offset_df_2[(offset_df_2['cred_bef']==1) & (offset_df_2['cred']==1)]
    fovs_network_draw_color_edge2(offset_df_2_cred,positions, maxweight,minweight,'um')
    maxweight,minweight,scale = fovs_network_draw_color_edge(offset_df_2_cred,positions,'nm')
    
    offset_df_2.to_csv(save_folder + '/offsets_1rnd.csv')
    identified_dup_spots_2.to_csv(save_folder + '/dup_spots_1rnd.csv')
    
    dup_ids = list(identified_dup_spots_2['id1'])+list(identified_dup_spots_2['id2'])
    spots_non_dup = spots.loc[~(spots['id'].isin(dup_ids))]
    
    
    positions['x']=positions['x'].astype(float)
    positions['y']=positions['y'].astype(float)
    minus_x = np.min(positions['x'])
    minus_y = np.min(positions['y'])
    if minus_x <0:
        positions['x']-=minus_x
        spots['global_x']-=minus_x*0.103

    if minus_y <0:
        positions['y']-=minus_y
        spots['global_y']-=minus_y*0.103
        
    spots.to_csv()
    positions.to_csv()    
    
    
    

