"""
    Copyright (C) 2018 The University of Sydney, Australia
    
    This program is free software; you can redistribute it and/or modify it under
    the terms of the GNU General Public License, version 2, as published by
    the Free Software Foundation.
    
    This program is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.
    
    You should have received a copy of the GNU General Public License along
    with this program; if not, write to Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

import pygplates
import numpy as np


def get_unique_plate_pairs_from_rotation_model(rotation_model,recon_time):
    # given a rotation model and a specifief reconstruction time, return 
    # a list of the unique plate pairs
    
    tree = rotation_model.get_reconstruction_tree(recon_time)
    edges = tree.get_edges()
    # Get a list of plate pairs
    tree_list = []
    for edge in edges:
        #if edge.get_parent_edge() is not None:
        #    tree_list.append((edge.get_fixed_plate_id(),edge.get_parent_edge().get_fixed_plate_id()))
        tree_list.append((edge.get_moving_plate_id(), edge.get_fixed_plate_id()))

    # get unique list (remove duplicates)
    uniq_plate_pairs_from_rotations = set(tree_list)

    return uniq_plate_pairs_from_rotations


def get_unique_plate_ids_from_reconstructed_features(reconstructed_features):
    # given a set of reconstructed features, return a list of the unique plate ids
    
    feature_plate_ids = []
    for reconstructed_feature in reconstructed_features:
        feature_plate_ids.append(reconstructed_feature.get_feature().get_reconstruction_plate_id())
    unique_plate_ids = set(feature_plate_ids)
    
    return unique_plate_ids


def GetPolygonCentroids(polygons):
    """
    Returns dict mapping plate IDs of static polygons to centroid of largest polygon
    associated with each plate ID.
    There is guaranteed to be a centroid for each plate ID.
    """
    
    is_resolved_topology = (type(polygons[0]) == pygplates.ResolvedTopologicalBoundary)
    
    centroid_dict = {}
    
    # Map plate ID to target area.
    target_polygon_area_dict = {}
    for polygon in polygons:
        plateid = polygon.get_feature().get_reconstruction_plate_id()
        if is_resolved_topology:
            area = polygon.get_resolved_geometry().get_area()
        else:
            area = polygon.get_reconstructed_geometry().get_area()
        
        # If first time adding area for plate ID or polygon area largest so far then
        # set new target area and also set new centroid.
        target_polygon_area = target_polygon_area_dict.get(plateid)
        if target_polygon_area is None or area > target_polygon_area:
            target_polygon_area_dict[plateid] = area
            if is_resolved_topology:
                centroid_dict[plateid] = ( # parentheses just so can put on next line...
                        polygon.get_resolved_geometry().get_boundary_centroid().to_lat_lon())
            else:
                centroid_dict[plateid] = ( # parentheses just so can put on next line...
                        polygon.get_reconstructed_geometry().get_boundary_centroid().to_lat_lon())
    
    return centroid_dict

################ UNUSED????? START
# function to get centroid from every polygon in the reconstructed static polygons
def GetPolygonCentroid(static_polygons,plateid):
    centroid = None
    target_polygon_area = 0
    for polygon in static_polygons:
        if polygon.get_feature().get_reconstruction_plate_id()==plateid:
            if polygon.get_reconstructed_geometry().get_area()>target_polygon_area:
                centroid = polygon.get_reconstructed_geometry().get_boundary_centroid().to_lat_lon()
                target_polygon_area = polygon.get_reconstructed_geometry().get_area()
    
    return centroid

# Alternatively, get centroids of topological polygons
def GetPlateCentroid(resolved_polygons,plateid):
    centroid = None
    for polygon in resolved_polygons:
        if polygon.get_feature().get_reconstruction_plate_id()==plateid:
            centroid = polygon.get_resolved_boundary().get_boundary_centroid().to_lat_lon()
    
    return centroid
################ UNUSED????? END


def get_root_static_polygon_plate_ids(reconstruction_tree, uniq_plates_from_static_polygons):
    """
    Find plates closest to root/anchor plate that are in the static polygons.
    The root/anchor plate (or even its children, and so on) may not be an actual static polygon.
    These are root static polygon plates - typically there's only one though.
    """
    
    # If anchor plate is in the static polygons then it's the only root plate.
    if reconstruction_tree.get_anchor_plate_id() in uniq_plates_from_static_polygons:
        return [reconstruction_tree.get_anchor_plate_id()]
    
    root_plates = []
    
    
    def traverse_sub_tree(edge):
        # If moving plate of current edge is in static polygons then
        # don't need to traverse child sub-tree.
        if edge.get_moving_plate_id() in uniq_plates_from_static_polygons:
            root_plates.append(edge.get_moving_plate_id())
            return
        
        # Recurse into the children sub-trees.
        for child_edge in edge.get_child_edges():
            traverse_sub_tree(child_edge)
    
    
    # Keep traversing towards the leaves until we find a moving plate ID that's in the static polygons.
    # Note that there can be multiple of these plate IDs.
    for anchor_plate_edge in reconstruction_tree.get_anchor_plate_edges():
        traverse_sub_tree(anchor_plate_edge)
    
    return root_plates


def patch_links_between_polygon(moving_plate, reconstruction_tree_edges_dict, uniq_plates_from_static_polygons):
    # For a given plate id, find the next highest plate id in the hierarchy
    # for which a geometry exists in the specified set of polygons,
    # otherwise return None.
    
    # Get edge in reconstruction tree associated with moving plate.
    edge_in_chain = reconstruction_tree_edges_dict.get(moving_plate)
    
    # If moving plate is anchor plate then there are no parent links.
    if edge_in_chain is None:
        return # Returns None
    
    # The first entry in the returned list is always the input moving plate.
    plate_chain_list = [moving_plate]
    
    # Traverse the reconstruction tree towards to the root/anchor plate.
    while True:
        plate_chain_list.append(edge_in_chain.get_fixed_plate_id())
        
        # Finished if fixed plate ID has a valid geometry.
        if edge_in_chain.get_fixed_plate_id() in uniq_plates_from_static_polygons:
            # We're guaranteed to get at least two list entries.
            return plate_chain_list
        
        # Move up the hierarachy towards the root/anchor plate.
        edge_in_chain = edge_in_chain.get_parent_edge()
        
        # If we've reached root/anchor plate but have not yet encountered
        # a fixed plate ID with a valid geometry then return None.
        if not edge_in_chain:
            break
    
    # By default a function returns None when reaches end.


def get_plate_chains(uniq_plates_from_static_polygons, reconstruction_tree):
    # given a list of plate ids found in static polygons, and a list of plate pairs
    # from a rotation tree, finds the linkages between plate geometries, and 
    # returns them as a list of lists.
    # where two plates (for which geometries exist) are linked by one or more
    # intermediate plate ids for which no geometry exists, the returned list will
    # contain three or more plate ids, where only the first and last entries
    # have valid geometries in the polygon set
    # where a plate is moving wrt to the spin axis, the list ends with zero
    
    # Map the moving plate ID of each edge (in reconstruction tree) to that edge for quick lookup.
    # Do it once here instead of inside each call to 'patch_links_between_polygon()'.
    reconstruction_tree_edges_dict = dict(
            (edge.get_moving_plate_id(), edge) # (key, value)
            for edge in reconstruction_tree.get_edges())
    
    chains = []
    for plate in uniq_plates_from_static_polygons: 
        # call function to patch links in plate chain
        chain = patch_links_between_polygon(
                plate,
                reconstruction_tree_edges_dict,
                uniq_plates_from_static_polygons)
        # If chain is None, means that no fixed plate geometry could be found for this moving plate.
        if chain:
            chains.append(chain)
        else:
            #print 'Root plate geometry with plate id %d' % plate
            continue
            
    return chains


def create_hierarchy_features(chains,reconstructed_polygons,tree_features=None,valid_time=None):
    #take plate chains and static polygons, and create a set of line features that 
    # join up the centroid points of polygons based on their linkage in the rotation
    # hierarchy. 
    # If tree_features in given as an existing list of features, the 
    # new features will be appended. Otherwise a new feature is created.
    # valid time (optional) can be given as a tuple
    
    if tree_features is None:
        tree_features = []

    polygon_centroids = GetPolygonCentroids(reconstructed_polygons)
    
    for chain in chains:

        # First and last plate IDs in a chain always correspond to a static polygon.
        # So p0 and p1 should not be None.
        p0 = polygon_centroids[chain[0]]
        p1 = polygon_centroids[chain[-1]]

        feature = pygplates.Feature()
        simple_line = pygplates.PolylineOnSphere([p0,p1])
        feature.set_geometry(simple_line.to_tessellated(np.radians(1)))
        feature.set_name(str(chain))
        if valid_time is not None:
            feature.set_valid_time(valid_time[0],valid_time[1])

        tree_features.append(feature)

    return tree_features
