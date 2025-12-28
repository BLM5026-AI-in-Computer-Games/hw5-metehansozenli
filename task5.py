#!/usr/bin/env python3

import argparse
import json
import random
import time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import List, Tuple

import folium
import networkx as nx
import numpy as np
import osmnx as ox
import pandas as pd
from ortools.constraint_solver import pywrapcp, routing_enums_pb2
from pyproj import Geod


@dataclass
class RegionCircle:
    region_id: int
    lat: float
    lon: float
    radius_m: float
    is_start: bool = False  # True for first region (no circle drawn)


@dataclass
class SolverRun:
    solver: str
    order: List[int]                           
    entry_points: List[Tuple[float, float]]    
    total_length_m: float
    runtime_s: float


_GEOD = Geod(ellps="WGS84")


def geodesic_distance(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    _, _, dist_m = _GEOD.inv(lon1, lat1, lon2, lat2)
    return float(dist_m)


def safe_shortest_path_length(G: nx.MultiDiGraph, u: int, v: int, weight: str = "length") -> float:
    if u == v:
        return 0.0
    try:
        return float(nx.shortest_path_length(G, u, v, weight=weight))
    except nx.NetworkXNoPath:
        return float("inf")


def safe_shortest_path(G: nx.MultiDiGraph, u: int, v: int, weight: str = "length") -> List[int]:
    if u == v:
        return [u]
    try:
        return nx.shortest_path(G, u, v, weight=weight)
    except nx.NetworkXNoPath:
        return []


def distinct_colors(n: int) -> List[str]:
    palette = [
        "#FF6B6B",  # Coral Red
        "#4ECDC4",  # Turquoise
        "#F7DC6F",  # Sunny Yellow
        "#7209B7",  # Purple
        "#52B788",  # Forest Green
        "#F8B739",  # Orange
        "#457B9D",  # Steel Blue
        "#E63946",  # Crimson
        "#06D6A0",  # Emerald
        "#F72585",  # Pink
        "#45B7D1",  # Sky Blue
        "#E76F51",  # Burnt Sienna
        "#98D8C8",  # Mint
        "#FB5607",  # Orange Red
        "#3A0CA3",  # Indigo
        "#FFA07A",  # Light Salmon
        "#2A9D8F",  # Teal
        "#BB8FCE",  # Lavender
        "#FFBE0B",  # Gold
        "#8338EC",  # Violet
        "#06FFA5",  # Neon Green
        "#A8DADC",  # Powder Blue
        "#F4A261",  # Sandy Brown
        "#4CC9F0",  # Bright Blue
        "#85C1E2",  # Baby Blue
    ]
    return [palette[i % len(palette)] for i in range(n)]


def generate_random_regions(place: str, n: int, radius_m: float, seed: int, G: nx.MultiDiGraph) -> List[RegionCircle]:
    random.seed(seed)
    np.random.seed(seed)

    gdf = ox.geocode_to_gdf(place)
    geom = gdf.geometry.iloc[0]

    # Collect candidate nodes inside the polygon
    candidates = []
    for node, data in G.nodes(data=True):
        lon = float(data.get("x"))
        lat = float(data.get("y"))
        pt = ox.utils_geo.Point(lon, lat)
        if geom.contains(pt):
            candidates.append(int(node))

    if len(candidates) < n:
        raise RuntimeError(f"Not enough road nodes inside '{place}' to sample {n} regions (found {len(candidates)}).")

    sampled = random.sample(candidates, n)
    regions = []
    for i, node in enumerate(sampled):
        lat = float(G.nodes[node]["y"])
        lon = float(G.nodes[node]["x"])
        is_start = (i == 0)  # First region is starting point
        regions.append(RegionCircle(region_id=i, lat=lat, lon=lon, radius_m=radius_m, is_start=is_start))
    return regions


# -----------------------------
# Circle intersection detection
# -----------------------------

def line_circle_intersection(p1_lat: float, p1_lon: float, p2_lat: float, p2_lon: float, center_lat: float, center_lon: float, radius_m: float) -> Tuple[float, float]:
    # Calculate distances from center
    dist_p1 = geodesic_distance(p1_lat, p1_lon, center_lat, center_lon)
    dist_p2 = geodesic_distance(p2_lat, p2_lon, center_lat, center_lon)
    
    # If no crossing, return closer point
    if (dist_p1 > radius_m and dist_p2 > radius_m) or (dist_p1 <= radius_m and dist_p2 <= radius_m):
        return (p1_lat, p1_lon) if abs(dist_p1 - radius_m) < abs(dist_p2 - radius_m) else (p2_lat, p2_lon)
    
    # Ensure p1 is outside, p2 is inside
    if dist_p1 <= radius_m:
        p1_lat, p2_lat = p2_lat, p1_lat
        p1_lon, p2_lon = p2_lon, p1_lon
    
    # Binary search for intersection point
    geod = Geod(ellps="WGS84")
    tolerance = 1.0  # 1 meter
    max_iterations = 50
    
    for _ in range(max_iterations):
        # Calculate midpoint along geodesic
        total_dist = geodesic_distance(p1_lat, p1_lon, p2_lat, p2_lon)
        if total_dist < 0.1:
            return (p1_lat, p1_lon)
        
        fwd_azimuth, _, _ = geod.inv(p1_lon, p1_lat, p2_lon, p2_lat)
        mid_lon, mid_lat, _ = geod.fwd(p1_lon, p1_lat, fwd_azimuth, total_dist / 2)
        
        dist_mid = geodesic_distance(mid_lat, mid_lon, center_lat, center_lon)
        
        if abs(dist_mid - radius_m) < tolerance:
            return (float(mid_lat), float(mid_lon))
        
        # Adjust search range
        if dist_mid > radius_m:
            p1_lat, p1_lon = mid_lat, mid_lon  # Mid is outside, move towards p2
        else:
            p2_lat, p2_lon = mid_lat, mid_lon  # Mid is inside, move towards p1
    
    # Return final approximation
    total_dist = geodesic_distance(p1_lat, p1_lon, p2_lat, p2_lon)
    fwd_azimuth, _, _ = geod.inv(p1_lon, p1_lat, p2_lon, p2_lat)
    final_lon, final_lat, _ = geod.fwd(p1_lon, p1_lat, fwd_azimuth, total_dist / 2)
    return (float(final_lat), float(final_lon))


def find_boundary_intersection(G: nx.MultiDiGraph, path_nodes: List[int], target_region: RegionCircle) -> Tuple[int, Tuple[float, float]]:
    if target_region.is_start:
        node = path_nodes[-1]
        node_lat = float(G.nodes[node]["y"])
        node_lon = float(G.nodes[node]["x"])
        return node, (node_lat, node_lon)
    
    # Walk through path and detect crossing
    for i in range(1, len(path_nodes)):
        prev_node = path_nodes[i - 1]
        curr_node = path_nodes[i]
        prev_lat = float(G.nodes[prev_node]["y"])
        prev_lon = float(G.nodes[prev_node]["x"])
        curr_lat = float(G.nodes[curr_node]["y"])
        curr_lon = float(G.nodes[curr_node]["x"])   
        prev_dist = geodesic_distance(prev_lat, prev_lon, target_region.lat, target_region.lon)
        curr_dist = geodesic_distance(curr_lat, curr_lon, target_region.lat, target_region.lon)
        
        # Crossing detected: prev outside, curr inside
        if prev_dist > target_region.radius_m and curr_dist <= target_region.radius_m:
            entry_point = line_circle_intersection(
                prev_lat, prev_lon,
                curr_lat, curr_lon,
                target_region.lat, target_region.lon,
                target_region.radius_m
            )
            return prev_node, entry_point
    
    # No crossing - find closest point to boundary
    best_node = path_nodes[-1]
    best_point = (float(G.nodes[best_node]["y"]), float(G.nodes[best_node]["x"]))
    min_dist_to_boundary = float('inf')
    
    for node in path_nodes:
        node_lat = float(G.nodes[node]["y"])
        node_lon = float(G.nodes[node]["x"])
        dist_to_center = geodesic_distance(node_lat, node_lon, target_region.lat, target_region.lon)
        dist_to_boundary = abs(dist_to_center - target_region.radius_m)
        
        if dist_to_boundary < min_dist_to_boundary:
            min_dist_to_boundary = dist_to_boundary
            best_node = node
            best_point = (node_lat, node_lon)
    
    return best_node, best_point


# -----------------------------
# Distance matrix for TSP ordering
# -----------------------------

def build_center_distance_matrix(G: nx.MultiDiGraph, regions: List[RegionCircle]) -> Tuple[np.ndarray, List[int]]:
    """Distance matrix between region centers"""
    center_nodes = [int(ox.distance.nearest_nodes(G, r.lon, r.lat)) for r in regions]
    n = len(regions)
    D = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(n):
            if i == j:
                D[i, j] = 0.0
            else:
                d = safe_shortest_path_length(G, center_nodes[i], center_nodes[j], weight="length")
                D[i, j] = d
    return D, center_nodes


class NearestNeighborSolver:
    """Greedy Nearest Neighbor TSP solver"""
    def solve(self, D: np.ndarray, start: int = 0) -> List[int]:
        n = D.shape[0]
        unvisited = set(range(n))
        tour = [start]
        unvisited.remove(start)
        while unvisited:
            cur = tour[-1]
            nxt = min(unvisited, key=lambda j: D[cur, j])
            tour.append(nxt)
            unvisited.remove(nxt)
        return tour


class ORToolsSolver:
    """OR-Tools TSP solver"""
    def __init__(self, time_limit_s: int = 10):
        self.time_limit_s = time_limit_s

    def solve(self, D: np.ndarray, start: int = 0) -> List[int]:
        n = D.shape[0]
        manager = pywrapcp.RoutingIndexManager(n, 1, start)
        routing = pywrapcp.RoutingModel(manager)

        scale = 1000.0

        def dist_cb(from_index: int, to_index: int) -> int:
            i = manager.IndexToNode(from_index)
            j = manager.IndexToNode(to_index)
            d = D[i, j]
            if not np.isfinite(d):
                return int(1e12)
            return int(d * scale)

        transit_idx = routing.RegisterTransitCallback(dist_cb)
        routing.SetArcCostEvaluatorOfAllVehicles(transit_idx)

        params = pywrapcp.DefaultRoutingSearchParameters()
        params.first_solution_strategy = routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
        params.local_search_metaheuristic = routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        params.time_limit.FromSeconds(int(self.time_limit_s))

        sol = routing.SolveWithParameters(params)
        if sol is None:
            return list(range(n))

        idx = routing.Start(0)
        order = []
        while not routing.IsEnd(idx):
            order.append(manager.IndexToNode(idx))
            idx = sol.Value(routing.NextVar(idx))
        return order


class GeneticAlgorithmSolver:
    """Genetic Algorithm TSP solver"""
    def __init__(
        self,
        population_size: int = 120,
        generations: int = 300,
        crossover_rate: float = 0.8,
        mutation_rate: float = 0.25,
        tournament_size: int = 5,
        elitism_count: int = 2,
        seed: int = 0,
    ):
        self.population_size = population_size
        self.generations = generations
        self.crossover_rate = crossover_rate
        self.mutation_rate = mutation_rate
        self.tournament_size = tournament_size
        self.elitism_count = elitism_count
        self.seed = seed

    def solve(self, D: np.ndarray, start: int = 0) -> List[int]:
        random.seed(self.seed)
        np.random.seed(self.seed)

        n = D.shape[0]
        cities = [i for i in range(n) if i != start]

        def fitness(ind: List[int]) -> float:
            tour = [start] + ind + [start]
            total = 0.0
            for a, b in zip(tour[:-1], tour[1:]):
                d = D[a, b]
                if not np.isfinite(d):
                    return float("inf")
                total += d
            return total

        def init_population() -> List[List[int]]:
            pop = []
            for _ in range(self.population_size):
                ind = cities.copy()
                random.shuffle(ind)
                pop.append(ind)
            return pop

        def tournament_select(pop: List[List[int]], scores: List[float]) -> List[int]:
            idxs = random.sample(range(len(pop)), self.tournament_size)
            best = min(idxs, key=lambda i: scores[i])
            return pop[best].copy()

        def order_crossover(p1: List[int], p2: List[int]) -> Tuple[List[int], List[int]]:
            size = len(p1)
            a, b = sorted(random.sample(range(size), 2))
            c1 = [None] * size
            c2 = [None] * size
            c1[a:b] = p1[a:b]
            c2[a:b] = p2[a:b]

            def fill(child, parent):
                s = set(x for x in child if x is not None)
                rem = [x for x in parent if x not in s]
                ri = 0
                for i in range(size):
                    if child[i] is None:
                        child[i] = rem[ri]
                        ri += 1
                return child

            return fill(c1, p2), fill(c2, p1)

        def swap_mutation(ind: List[int]) -> List[int]:
            if len(ind) < 2:
                return ind
            i, j = random.sample(range(len(ind)), 2)
            ind2 = ind.copy()
            ind2[i], ind2[j] = ind2[j], ind2[i]
            return ind2

        pop = init_population()

        for _ in range(self.generations):
            scores = [fitness(ind) for ind in pop]
            elite_idxs = np.argsort(scores)[: self.elitism_count]
            new_pop = [pop[i].copy() for i in elite_idxs]

            while len(new_pop) < self.population_size:
                p1 = tournament_select(pop, scores)
                p2 = tournament_select(pop, scores)

                if random.random() < self.crossover_rate and len(p1) >= 2:
                    c1, c2 = order_crossover(p1, p2)
                else:
                    c1, c2 = p1, p2

                if random.random() < self.mutation_rate:
                    c1 = swap_mutation(c1)
                if random.random() < self.mutation_rate:
                    c2 = swap_mutation(c2)

                new_pop.extend([c1, c2])

            pop = new_pop[: self.population_size]

        scores = [fitness(ind) for ind in pop]
        best = pop[int(np.argmin(scores))]
        return [start] + best


def compute_tour( G: nx.MultiDiGraph, regions: List[RegionCircle], order: List[int], center_nodes: List[int]) -> Tuple[List[Tuple[float, float]], List[List[int]], float]:
    entry_points = []
    leg_paths = []
    total_length = 0.0
    start_node = center_nodes[order[0]]
    start_lat = float(G.nodes[start_node]["y"])
    start_lon = float(G.nodes[start_node]["x"])
    current_entry_point = (start_lat, start_lon)
    
    entry_points.append(current_entry_point)

    for i in range(len(order)):
        next_idx = (i + 1) % len(order)
        next_region_id = order[next_idx]
        next_region = regions[next_region_id]
        next_center = center_nodes[next_region_id]

        start_node = int(ox.distance.nearest_nodes(
            G, 
            current_entry_point[1],  # lon
            current_entry_point[0]   # lat
        ))
        
        # Compute path from start_node to next region's center
        path_to_center = safe_shortest_path(G, start_node, next_center, weight="length")
        
        if not path_to_center:
            raise RuntimeError(f"No path from node {start_node} to region {next_region_id} center")
        
        if path_to_center[0] != start_node:
            path_to_center = [start_node] + path_to_center
        
        entry_node, next_entry_point = find_boundary_intersection(G, path_to_center, next_region)
        
        # Truncate path to entry node
        try:
            entry_idx = path_to_center.index(entry_node)
            actual_path = path_to_center[:entry_idx + 1]
        except ValueError:
            actual_path = path_to_center
            entry_node = path_to_center[-1]
        
        if len(actual_path) < 2:
            actual_path = [start_node, entry_node] if start_node != entry_node else [start_node]
        
        # Calculate leg length from current_entry_point to next_entry_point:
        # 1. Distance from current_entry_point to start_node
        start_node_lat = float(G.nodes[start_node]["y"])
        start_node_lon = float(G.nodes[start_node]["x"])
        dist_from_entry_to_start = geodesic_distance(
            current_entry_point[0], current_entry_point[1],
            start_node_lat, start_node_lon
        )
        
        # 2. Path length along road network
        if len(actual_path) >= 2:
            path_length = float(nx.path_weight(G, actual_path, weight="length"))
        else:
            path_length = 0.0
        
        # 3. Distance from entry_node to next_entry_point
        entry_node_lat = float(G.nodes[entry_node]["y"])
        entry_node_lon = float(G.nodes[entry_node]["x"])
        dist_from_end_to_entry = geodesic_distance(
            entry_node_lat, entry_node_lon,
            next_entry_point[0], next_entry_point[1]
        )
        
        # Total leg length
        leg_length = dist_from_entry_to_start + path_length + dist_from_end_to_entry
        total_length += leg_length
        
        leg_paths.append(actual_path)
        entry_points.append(next_entry_point)
        current_entry_point = next_entry_point
    
    entry_points = entry_points[:-1]
    
    return entry_points, leg_paths, total_length


# -----------------------------
# Route geometry extraction
# -----------------------------

def extract_route_geometry(G: nx.MultiDiGraph, path_nodes: List[int]) -> List[Tuple[float, float]]:
    """Extract detailed route geometry (lat, lon) from path nodes using edge geometries"""
    if len(path_nodes) < 2:
        # Single node, return its coordinates
        if len(path_nodes) == 1:
            node = path_nodes[0]
            return [(float(G.nodes[node]["y"]), float(G.nodes[node]["x"]))]
        return []
    
    coords = []
    
    for i in range(len(path_nodes) - 1):
        u = path_nodes[i]
        v = path_nodes[i + 1]
        
        # Add starting node coordinates
        if i == 0:
            coords.append((float(G.nodes[u]["y"]), float(G.nodes[u]["x"])))
        
        # Check if edge has geometry attribute (from OSMnx simplification)
        edge_data = None
        if G.has_edge(u, v):
            # Get the first edge (key=0) if it exists
            edges = G[u][v]
            if edges:
                edge_data = edges[list(edges.keys())[0]]
        
        if edge_data and 'geometry' in edge_data:
            # Use the detailed geometry
            geom = edge_data['geometry']
            for point in list(geom.coords)[1:]:
                coords.append((point[1], point[0])) 
        else:
            # No geometry, just add the end node
            coords.append((float(G.nodes[v]["y"]), float(G.nodes[v]["x"])))
    
    return coords


# -----------------------------
# Folium visualization
# -----------------------------

def visualize_tour(
    place: str,
    G: nx.MultiDiGraph,
    regions: List[RegionCircle],
    run: SolverRun,
    leg_paths: List[List[int]],
    out_html: Path,
):
    """Create interactive folium map"""
    center_lat = float(regions[run.order[0]].lat)
    center_lon = float(regions[run.order[0]].lon)

    fmap = folium.Map(location=[center_lat, center_lon], zoom_start=12, tiles="cartodbpositron")

    # Draw regions
    for r in regions:
        if not r.is_start:
            # Draw circle for non-start regions with fill
            folium.Circle(
                location=[r.lat, r.lon],
                radius=r.radius_m,
                color="#536476",
                weight=2,
                fill=True,
                fill_color="#3498db",
                fill_opacity=0.35,
                opacity=0.8,
                tooltip=f"Region {r.region_id} (r={int(r.radius_m)}m)",
            ).add_to(fmap)
        
        # Draw center marker
        icon_html = f"""<div style="font-size: 12pt; color: #2C3E50; font-weight: 600;"><b>{r.region_id}</b></div>"""
        if r.is_start:
            icon_html = f"""<div style="font-size: 12pt; color: #E74C3C; font-weight: 700;"><b>START {r.region_id}</b></div>"""
        
        folium.Marker(
            location=[r.lat, r.lon],
            icon=folium.DivIcon(html=icon_html),
        ).add_to(fmap)

    # Mark exact entry points on circle boundaries
    for idx_in_order, entry_point in enumerate(run.entry_points):
        lat, lon = entry_point
        region_id = run.order[idx_in_order]
        region = regions[region_id]
        
        # Calculate distance from entry point to center (should be ~radius for non-start regions)
        dist_to_center = geodesic_distance(lat, lon, region.lat, region.lon)
        dist_to_boundary = abs(dist_to_center - region.radius_m)
        
        color = "#B42727" if idx_in_order == 0 else "#60CE60"
        folium.CircleMarker(
            location=[lat, lon],
            radius=8,
            color=color,
            weight=3,
            fill=True,
            fill_color="white",
            fill_opacity=1.0,
            tooltip=f"EXACT Entry for Region {region_id}<br/>Distance to center: {dist_to_center:.2f}m<br/>Distance to boundary: {dist_to_boundary:.2f}m<br/>Radius: {region.radius_m:.2f}m",
        ).add_to(fmap)

    # Draw route legs with proper geometry - each leg starts from previous entry_point and ends at next entry_point
    colors = distinct_colors(len(leg_paths))
    for i, path in enumerate(leg_paths):
        # Extract detailed route geometry from the node path
        coords = extract_route_geometry(G, path)

        prev_entry_idx = i  # Current leg goes FROM entry i TO entry i+1
        next_entry_idx = (i + 1) % len(run.entry_points)
        
        # Build the complete leg coordinates:
        leg_coords = [run.entry_points[prev_entry_idx]]
        
        # Add all the intermediate road geometry
        for coord in coords:
            leg_coords.append(coord)
        
        # End with next entry_point (exact)
        leg_coords.append(run.entry_points[next_entry_idx])
        
        if leg_coords and len(leg_coords) >= 2:
            folium.PolyLine(leg_coords, color=colors[i], weight=5, opacity=0.8, 
                          tooltip=f"Leg {i}: Region {run.order[prev_entry_idx]} → {run.order[next_entry_idx]}").add_to(fmap)

    html = f"""
    <div style="
      position: fixed;
      bottom: 20px; left: 20px;
      z-index: 9999;
      background: white;
      padding: 12px 14px;
      border: 1px solid #999;
      border-radius: 8px;
      box-shadow: 2px 2px 8px rgba(0,0,0,0.2);
      font-size: 12px;
      max-width: 360px;">
      <b>Task 5 — {run.solver}</b><br/>
      Place: {place}<br/>
      Regions: {len(regions)}<br/>
      Total length: <b>{run.total_length_m/1000:.3f} km</b><br/>
      Runtime: {run.runtime_s:.3f} s<br/>
    </div>
    """
    fmap.get_root().html.add_child(folium.Element(html))
    fmap.save(str(out_html))


def run_solver_pipeline(
    solver_name: str,
    solver,
    G: nx.MultiDiGraph,
    regions: List[RegionCircle],
    D: np.ndarray,
    center_nodes: List[int],
    start: int = 0
) -> Tuple[SolverRun, List[List[int]]]:
    t0 = time.perf_counter()
    order = solver.solve(D, start=start)
    
    entry_points, leg_paths, total_length = compute_tour(
        G, regions, order, center_nodes
    )
    
    t1 = time.perf_counter()
    runtime = t1 - t0
    
    run = SolverRun(
        solver=solver_name,
        order=order,
        entry_points=entry_points,
        total_length_m=total_length,
        runtime_s=runtime,
    )
    return run, leg_paths


def main():
    parser = argparse.ArgumentParser(description="Task 5: TSPN with Circle Entry Detection")
    parser.add_argument("--place", type=str, default="Nilüfer, Bursa, Türkiye")
    parser.add_argument("--regions", type=int, default=10)
    parser.add_argument("--radius-m", type=float, default=250.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--graph-cache", type=str, default="results/nilufer_drive.graphml")
    parser.add_argument("--or-time-limit", type=int, default=12)
    parser.add_argument("--ga-pop", type=int, default=120)
    parser.add_argument("--ga-gen", type=int, default=300)
    parser.add_argument("--outdir", type=str, default="results")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\nTask 5: TSP with Neighborhoods")
    print(f"Place: {args.place}")
    print(f"Regions: {args.regions}, Radius: {args.radius_m}m\n")

    # Load or download graph
    graph_cache = Path(args.graph_cache)
    graph_cache.parent.mkdir(parents=True, exist_ok=True)

    if graph_cache.exists():
        print("Loading cached graph...")
        G = ox.load_graphml(graph_cache)
    else:
        print("Downloading road network...")
        G = ox.graph_from_place(args.place, network_type="drive", simplify=True)
        ox.save_graphml(G, graph_cache)

    if not any("length" in data for _, _, data in G.edges(data=True)):
        G = ox.add_edge_lengths(G)

    # Generate regions
    print("Generating regions...")
    regions = generate_random_regions(args.place, args.regions, args.radius_m, args.seed, G)

    if len(regions) < 3:
        raise ValueError("Need at least 3 regions")

    # Build distance matrix
    print("Building distance matrix...")
    D, center_nodes = build_center_distance_matrix(G, regions)

    # Run three solvers
    print("\nRunning solvers...")
    solvers = [
        ("Nearest Neighbor", NearestNeighborSolver()),
        ("OR-Tools", ORToolsSolver(time_limit_s=args.or_time_limit)),
        ("Genetic Algorithm", GeneticAlgorithmSolver(
            population_size=args.ga_pop,
            generations=args.ga_gen,
            seed=args.seed
        )),
    ]

    results = []
    for solver_name, solver in solvers:
        print(f"  Running {solver_name}...")
        run, leg_paths = run_solver_pipeline(
            solver_name, solver, G, regions, D, center_nodes, start=args.start
        )
        
        # Save results
        slug = solver_name.split()[0].lower()
        
        # JSON
        json_path = outdir / f"paths_{slug}.json"
        with open(json_path, "w", encoding="utf-8") as f:
            json.dump({
                "solver": run.solver,
                "place": args.place,
                "regions": [asdict(r) for r in regions],
                "order": run.order,
                "entry_points": run.entry_points,
                "total_length_m": run.total_length_m,
                "runtime_s": run.runtime_s,
                "leg_paths_nodes": leg_paths,
            }, f, ensure_ascii=False, indent=2)
        
        # HTML map
        html_path = outdir / f"task5_{slug}.html"
        visualize_tour(args.place, G, regions, run, leg_paths, html_path)
        
        results.append({
            "solver": run.solver,
            "total_length_km": run.total_length_m / 1000.0,
            "runtime_s": run.runtime_s,
            "order": run.order,
            "html_map": str(html_path),
        })
        
        print(f"    ✓ Length: {run.total_length_m/1000:.3f} km, Time: {run.runtime_s:.3f}s")

    # Save summary
    df = pd.DataFrame(results)
    df.to_csv(outdir / "summary.csv", index=False)

    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(df[["solver", "total_length_km", "runtime_s"]].to_string(index=False))
    print(f"\nOutputs saved to: {outdir.resolve()}")
    print("=" * 70)


if __name__ == "__main__":
    main()
