import math
import pandas as pd
import numpy as np
import osmnx as ox
import networkx as nx
import time


def haversine(theta):
# Υπολογισμός του θ για το Great Circle

    return math.sin(theta / 2) ** 2

def great_circle_distance(lat1, lon1, lat2, lon2):
# Υπολογισμός της απόστασης Great Circle

    r = 6371  # Ακτίνα της Γης σε χιλιόμετρα
    lat1_rad = math.radians(lat1)
    lon1_rad = math.radians(lon1)
    lat2_rad = math.radians(lat2)
    lon2_rad = math.radians(lon2)

    delta_lat = lat2_rad - lat1_rad
    delta_lon = lon2_rad - lon1_rad

    a = haversine(delta_lat) + math.cos(lat1_rad) * math.cos(lat2_rad) * haversine(delta_lon)
    c = 2 * r * math.asin(math.sqrt(a))*1000

    return c

def find_projection_point(point_lat, point_lon, road_segment):
# Βρίσκει την προβολή (συντεταγμένες) του στίγματος GPS επί της ευθείας του οδικού τμήματος

    lat1, lon1 = road_segment[0]
    lat2, lon2 = road_segment[1]

    d_la2_la1 = lat2 - lat1
    d_lo2_lo1 = lon2 - lon1

    d_pla_la1 = point_lat - lat1
    d_plo_lo1 = point_lon - lon1

    coef = ((d_plo_lo1 * d_lo2_lo1 + d_pla_la1 * d_la2_la1) / (d_lo2_lo1 ** 2 + d_la2_la1 ** 2))

    proj_lat = lat1 + coef * (lat2 - lat1)
    proj_lon = lon1 + coef * (lon2 - lon1)
    projection = (proj_lat, proj_lon)

    return projection

def is_projection_outside_segment(projected_lat, projected_lon, road_segment):
# Ελέγχει αν η προβολή του GPS στίγματος είναι εντός του οδικού τμήματος

    lat1, lon1 = road_segment[0]
    lat2, lon2 = road_segment[1]

    min_lat = min(lat1, lat2)
    max_lat = max(lat1, lat2)
    min_lon = min(lon1, lon2)
    max_lon = max(lon1, lon2)

    if projected_lat < min_lat or projected_lat > max_lat or \
            projected_lon < min_lon or projected_lon > max_lon :
       return True
    else:
       return False

def select_point_of_segment(lat_point, lon_point, road_segment):
# Η προβολή είναι εκτός του οδικού τμήματος και επιλέγει το κοντινότερο άκρο του οδικού τμήματος

    lat1, lon1 = road_segment[0]
    lat2, lon2 = road_segment[1]

    # Υπολογισμός απόστασης από την αρχή του τμήματος δρόμου
    distance_start = great_circle_distance(lat1, lon1, lat_point, lon_point)

    # Υπολογισμός απόστασης από το τέλος του τμήματος δρόμου
    distance_end = great_circle_distance(lat2, lon2, lat_point, lon_point)

    # Επιλογή της νέας θέσης προβολής ως το άκρο με τη μικρότερη απόσταση
    if distance_start < distance_end:
        projected_lat, projected_lon = lat1, lon1
    else:
        projected_lat, projected_lon = lat2, lon2

    projection = (projected_lat, projected_lon)

    return projection

def shortest_path (place, start_point, end_point):
# Βρίσκει τη απόσταση (συντομότερη διαδρομή) μεταξύ δύο σημείων στο χάρτη

    switch = 1  # Με switch = 1 η απόσταση υπολογίζεται μέσω great circle, με switch = 2 μέσω OSMnx
    start_lat = start_point[0]
    start_lon = start_point[1]
    end_lat = end_point[0]
    end_lon = end_point[1]

# Μέθοδος great circle - υπολογισμός απόστασης μεταξύ 2 σημείων
    if switch == 1 :
        total_distance = great_circle_distance(start_lat, start_lon, end_lat, end_lon)

    else :
# Μέθοδος OSMNX - εύρεση της συντομότερης οδικής διαδρομής και υπολογισμός της απόστασης

        G = ox.graph_from_place(place, network_type='drive', which_result=1)

# Βρίσκουμε τους κοντινότερους κόμβους στα δίκτυα δρόμου
        start_node = ox.distance.nearest_nodes(G, start_lon, start_lat)
        end_node = ox.distance.nearest_nodes(G, end_lon, end_lat)

# Υπολογίζουμε τη συντομότερη διαδρομή με βάση τους κόμβους εκκίνησης και προορισμού
        shortest_path = nx.shortest_path(G, start_node, end_node, weight='length')

# Αποθηκεύουμε τους αναγνωριστικούς αριθμούς των κόμβων στη λίστα node_ids
        node_ids = shortest_path

# Μετατρέπουμε το γράφημα σε GeoDataFrames
        nodes = ox.graph_to_gdfs(G, nodes=True, edges=False)

# Επιλέγουμε τους κόμβους της συντομότερης διαδρομής
        selected_nodes = nodes[nodes.index.isin(node_ids)]

# Υπολογίζουμε τη συνολική απόσταση της συντομότερης διαδρομής
        total_distance = 0

# Διατρέχουμε τους κόμβους της συντομότερης διαδρομής
        for i in range(len(shortest_path) - 1):
    # Παίρνουμε το ακμή μεταξύ διαδοχικών κόμβων
            edge = G[shortest_path[i]][shortest_path[i+1]]
    # Προσθέτουμε το μήκος της ακμής στη συνολική απόσταση
            total_distance += edge[0]['length']

    return total_distance

def viterbi_map_matching(num_observations, candidate_li, emit_li, trans_li):
    """
    Υλοποιεί τον αλγόριθμο Viterbi για Αντιστοίχιση Χάρτη (Map Matching).

    Args:
    - num_observations: Το πλήθος των στιγμάτων GPS για τη συγκεκριμένη διαδρομή.
    - candidate_li: H λίστα με τις επιμέρους λίστες των υποψήφιων οδικών τμημάτων ανά στίγμα GPS για τη συγκεκριμένη
      διαδρομή.
    - emit_li: H λίστα με τις επιμέρους λίστες των πιθανοτήτων εκπομπής για τα υποψήφια οδικά τμήματα ανά στίγμα GPS για
      τη συγκεκριμένη διαδρομή.
    - trans_li: H λίστα με τις επιμέρους λίστες των πιθανοτήτων μετάβασης για τα υποψήφια οδικά τμήματα ανά στίγμα GPS
      (σε σχέση με τα υποψήφια οδικά τμήματα του προηγούμενου στίγματος) για τη συγκεκριμένη διαδρομή.

    Returns:
    - Τη λίστα path που περιέχει τις αντιστοιχίσεις των οδικών τμημάτων (τα ID τους) στα δοθέντα GPS στίγματα.
    """

    # Αρχικοποίηση του πίνακα Viterbi

    viterbi_prob = []
    path_list = []

    # Αρχικοποίηση για t = 0

    emit_prob = emit_li[0].split()
    candidate_items = candidate_li[0].split()
    num_cand = len(emit_prob)
    for cand in range(num_cand):
        viterbi_prob.append(float(emit_prob[cand]))
        path_list.append(candidate_items[cand])
    prev_cand = num_cand

    # Αναδρομικός υπολογισμός

    for t in range(1, num_observations):
        emit_prob = emit_li[t].split()
        candidate_items = candidate_li[t].split()
        num_cand = len(emit_prob)
        prev_viterbi_prob = viterbi_prob
        viterbi_prob = []
        prev_path_list = path_list
        path_list = []
        trans = 0
        for cand in range(num_cand):
            max_prob = 0
            max_backpointer = -1
            emission_prob = float(emit_prob[cand])
            trans_prob = trans_li[t - 1].split()
            for pr_cand in range(0, prev_cand):
                transition_prob = float(trans_prob[trans])
                trans = trans + 1
                viterbi_probability = prev_viterbi_prob[pr_cand]
                prob = viterbi_probability * transition_prob * emission_prob
                if prob > max_prob:
                    max_prob = prob
                    max_backpointer = pr_cand

            viterbi_prob.append(max_prob)
            cand_item = prev_path_list[max_backpointer] + ' ' + candidate_items[cand]
            path_list.append(cand_item)
        prev_cand = num_cand

        # Κανονικοποίηση των πιθανοτήτων viterbi_prob στο 1

        nsum = 0
        old_viterbi_prob = viterbi_prob
        viterbi_prob = []
        for cand in range(num_cand):
            nsum = nsum + old_viterbi_prob[cand]
        for cand in range(num_cand):
            viterbi_prob.append(old_viterbi_prob[cand] / nsum)

    # Βρίσκουμε την τελική κατάσταση με τη μεγαλύτερη πιθανότητα και το αντίστοιχο path

    final_state = np.argmax(viterbi_prob)
    path = path_list[final_state]

    return path


#---------------------------------------------------------------------------
"""
# ΤΟ ΚΕΝΤΡΙΚΟ ΠΡΟΓΡΑΜΜΑ 
Φάση 1. Υπολογισμός πιθανοτήτων αρχικής κατάστασης και πιθανοτήτων εκπομπής

"""
start_time = time.time()


excel_data_routes = pd.read_excel('Riga_Routes.xlsx')
count_routes = len(excel_data_routes)

excel_data_gps = pd.read_excel('GPS_Data.xlsx', dtype=object)
count_gps = len(excel_data_gps)
output_name = excel_data_gps['SelNode'].iloc[0]

sigma = 4.07
del_list = []
counter = 0


# Επιλέγονται τα υποψήφια οδικά τμήματα ανά GPS στίγμα και υπολογίζονται οι πιθανότητες

for row in range(count_gps) :
    candidates_id = []      # το ID του υποψήφιου οδικού τμήματος για το συγκεκριμένο στίγμα
    candidates_proj = []    # οι συντεταγμένες της προβολής του στίγματος στο υποψήφιο οδικό τμήμα
    candidates_pgc = []     # η GC απόσταση |Zt-Xt,i|
    candidates_emit = []    # οι πιθανότητες εκπομπής για κάθε υποψήφιο
    sum_emit_prob = 0

    pre_cand_id = excel_data_gps.iloc[row, 7].split()
    count_cand = len(pre_cand_id)
    point_lat = excel_data_gps.iloc[row, 3]
    point_lon = excel_data_gps.iloc[row, 4]

    for j in range(count_cand) :
        route = int(pre_cand_id[j])
        geometry = excel_data_routes.iloc[route-1, 4].split()
        count_segments = len(geometry)
        i = 0
        gps_segm_dist = 1000
        while i < (count_segments - 3) :
            road_segment = [(float(geometry[i]), float(geometry[i + 1])), (float(geometry[i + 2]), float(geometry[i + 3]))]

            projection_point = find_projection_point(point_lat, point_lon, road_segment)
            if is_projection_outside_segment(projection_point[0], projection_point[1], road_segment):
                projection_point = select_point_of_segment(point_lat, point_lon, road_segment)

            segm_dist = great_circle_distance(point_lat, point_lon, projection_point[0], projection_point[1])
            if segm_dist < gps_segm_dist :
               gps_segm_dist = segm_dist
               gps_proj_coord = projection_point
            i = i + 2

        if gps_segm_dist < 200:
            candidates_id.append(route)
            candidates_proj.append(gps_proj_coord[0])
            candidates_proj.append(gps_proj_coord[1])
            candidates_pgc.append(gps_segm_dist)
            emit_prob = (math.exp(1) ** (-0.5 * (gps_segm_dist / sigma) ** 2)) / (np.sqrt(math.pi * 2) * sigma)
            sum_emit_prob = sum_emit_prob + emit_prob
            candidates_emit.append(emit_prob)


    # Διαδικασία normalization πιθανοτήτων εκπομπής, αφαιρούνται πιθανότητες <= 0,0001

    norm_cand_emit = candidates_emit
    count_emit_prob = len(norm_cand_emit)
    norm_route = candidates_id
    norm_gps_proj = candidates_proj
    norm_pgc = candidates_pgc

    candidates_id = ''
    candidates_proj = ''
    candidates_pgc = ''
    candidates_emit = ''

    for i in range(count_emit_prob) :
        if sum_emit_prob > 0 :
            norm_emit_prob = norm_cand_emit[i] / sum_emit_prob
            if norm_emit_prob > 0.0001 :

                candidates_emit = candidates_emit + str(norm_emit_prob) + ' '
                candidates_id = candidates_id + str(norm_route[i]) + ' '
                candidates_proj = candidates_proj + str(norm_gps_proj[2*i]) + ' '
                candidates_proj = candidates_proj + str(norm_gps_proj[2 * i+1]) + ' '
                candidates_pgc = candidates_pgc + str(norm_pgc[i]) + ' '

    excel_data_gps.iloc[row, 8] = candidates_id
    excel_data_gps.iloc[row, 9] = candidates_proj
    excel_data_gps.iloc[row, 10] = candidates_pgc
    excel_data_gps.iloc[row, 11] = candidates_emit

    if not candidates_id:
        del_list.append(row)

    counter = counter + 1
    if counter % 10000 == 0 :
        print('Αθροιστής :', counter)


# Διαγραφή GPS στιγμάτων χωρίς υποψήφια οδικά τμήματα

if del_list :
    excel_data_gps.drop(del_list, inplace=True)
    count_gps = len(excel_data_gps)
    print('GPS στίγματα χωρίς υποψήφιους δρόμους: ', len(del_list))

    # Υπολογισμός great circle απόστασης 2 διαδοχικών GPS στιγμάτων |Zt-Zt+1|
    for row in range(count_gps):
        if row+1 < count_gps and (excel_data_gps.iloc[row, 2] == excel_data_gps.iloc[row+1, 2]):
            point_lat1 = excel_data_gps.iloc[row, 3]
            point_lon1 = excel_data_gps.iloc[row, 4]
            point_lat2 = excel_data_gps.iloc[row+1, 3]
            point_lon2 = excel_data_gps.iloc[row+1, 4]
            excel_data_gps.iloc[row, 6] = great_circle_distance(point_lat1, point_lon1, point_lat2, point_lon2)
        else :
            excel_data_gps.iloc[row, 6] = ''

elapsed_time1 = time.time() - start_time
print("Η Φάση 1 εκτελέστηκε σε", elapsed_time1, "δευτερόλεπτα.")


start_time2 = time.time()


"""
Φάση 2. Υπολογισμός των πιθανοτήτων μετάβασης

"""

#    Ορίζουμε την τοποθεσία "Riga, Latvia"
place ="Riga, Latvia"

for row in range(count_gps) :
    candidates_dist = ''  # η απόσταση |Xt,i-Xt+1,j|
    dt_list = []
    if row+1 < count_gps :
        vechicle_id = excel_data_gps.iloc[row, 2]
        vechicle_id_next = excel_data_gps.iloc[row+1, 2]
        if vechicle_id == vechicle_id_next :
            gps_dist = excel_data_gps.iloc[row, 6]
            candidates_proj = excel_data_gps.iloc[row, 9].split()
            count_proj = int(len(candidates_proj) / 2)
            candidates_proj_next = excel_data_gps.iloc[row+1, 9].split()
            count_proj_next = int(len(candidates_proj_next) / 2)
            for i in range(count_proj) :
                lat_curr = float(candidates_proj[i*2])
                lon_curr = float(candidates_proj[i*2+1])
                start_cand = (lat_curr, lon_curr)
                for j in range(count_proj_next) :
                    lat_next = float(candidates_proj_next[j*2])
                    lon_next = float(candidates_proj_next[j*2+1])
                    end_cand = (lat_next, lon_next)
                    cand_dist = shortest_path(place, start_cand, end_cand)
                    dt_list.append(np.abs(gps_dist-cand_dist))
            beta = np.median(dt_list) / (math.log(2))

            if beta == 0 :
                beta = 0.000001

            trans = 0
            sum_trans_prob = 0
            cand_trans_prob = []
            for i in range(count_proj) :
                for j in range(count_proj_next) :
                    trans_prob = (math.exp(1) ** (-dt_list[trans] / beta)) / beta
                    cand_trans_prob.append(trans_prob)
                    sum_trans_prob = sum_trans_prob + trans_prob
                    trans = trans + 1

            # Διαδικασία normalization πιθανοτήτων μετάβασης
            trans = 0
            candidates_trans_prob = ''
            for i in range(count_proj):
                for j in range(count_proj_next):
                    norm_trans_prob = cand_trans_prob[trans] / sum_trans_prob
                    candidates_trans_prob = candidates_trans_prob + str(norm_trans_prob) + ' '
                    trans = trans + 1

            excel_data_gps.iloc[row, 12] = candidates_trans_prob


elapsed_time2 = time.time() - start_time2
print("Η Φάση 2 εκτελέστηκε σε", elapsed_time2, "δευτερόλεπτα.")

"""
# Φάση 3. Επιλέγονται τα οδικά τμήματα και τα διορθωμένα στίγματα στο χάρτη. Αλγόριθμος Viterbi

"""

start_time3 = time.time()

# Μετράει ανά όχημα (διαδρομή) τις παρατηρήσεις (στίγματα GPS) και τα φυλάει στη λίστα obs_vechicle

i = 0
count_ov_start = 0
count_ov_end = 0
list_ov_start = []
list_ov_end = []

while i < count_gps:
    vech = excel_data_gps.iloc[i, 2]
    list_ov_start.append(count_ov_start)
    count_ov_end = count_ov_start
    while i < count_gps and vech == excel_data_gps.iloc[i, 2] :
        count_ov_end = count_ov_end + 1
        i = i + 1
    list_ov_end.append(count_ov_end)
    count_ov_start = count_ov_end

no_vech = len(list_ov_end)  # πλήθος οχημάτων (και διαδρομών)

# Υπολογίζει ανά όχημα την συντομότερη διαδρομή (Viterbi) με βάση τα δεδομένα στίγματα GPS

for i in range(no_vech):

    excel_stili = 8
    candidate_list = excel_data_gps.iloc[list_ov_start[i]:list_ov_end[i], excel_stili].tolist()
    excel_stili = 9
    projection_list = excel_data_gps.iloc[list_ov_start[i]:list_ov_end[i], excel_stili].tolist()
    excel_stili = 11
    emit_list = excel_data_gps.iloc[list_ov_start[i]:list_ov_end[i], excel_stili].tolist()
    excel_stili = 12
    trans_list = excel_data_gps.iloc[list_ov_start[i]:list_ov_end[i], excel_stili].tolist()
    no_obs_vechicle = ((list_ov_end[i]) - list_ov_start[i])

    path = viterbi_map_matching(no_obs_vechicle, candidate_list, emit_list, trans_list)

    # Μεταφέρει τα ID των επιλεγμένων οδικών τμημάτων και τις συντεταγμένες των επιλεγμένων προβολών στο excel αρχείο

    excel_path = path.split()
    no_excel_path = 0
    for j in range(list_ov_start[i], list_ov_end[i]):
        variable1 = excel_path[no_excel_path]  # το ID του επιλεγμένου οδικού τμήματος για το συγκεκριμένο στίγμα
        candidate_li = candidate_list[no_excel_path].split()
        cand_index = candidate_li.index(variable1)
        proj_li = projection_list[no_excel_path].split()
        variable2 = proj_li[cand_index * 2] + ' ' + proj_li[cand_index * 2 + 1]  # οι συντεταγμένες της προβολής του στίγματος
        excel_data_gps.iloc[j, 13] = variable1
        excel_data_gps.iloc[j, 14] = variable2

        no_excel_path = no_excel_path + 1

# Καταχώρηση στο excel αρχείο ΧΧΧΧ_Results

excel_data_gps.to_excel(f'{output_name}_Results.xlsx', index=False)

elapsed_time3 = time.time() - start_time2
print("Η Φάση 3 εκτελέστηκε σε", elapsed_time3, "δευτερόλεπτα.")

print("Το συνολικό πρόγραμμα εκτελέστηκε σε", elapsed_time1 + elapsed_time2 + elapsed_time3, "δευτερόλεπτα.")

