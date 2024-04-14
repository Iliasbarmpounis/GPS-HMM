import osmnx as ox
import pandas as pd
from rtree import index
import math
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
    delta_lat = abs(lat2_rad - lat1_rad)
    delta_lon = abs(lon2_rad - lon1_rad)
    a = haversine(delta_lat) + math.cos(lat1_rad) * math.cos(lat2_rad) * haversine(delta_lon)
    c = 2 * r * math.asin(math.sqrt(a))*1000

    return c

def rearrange_coordinates(line_string):
# Μορφοποίηση της στήλης geometry

    if 'LINESTRING' in str(line_string):
        pairs = line_string.replace("LINESTRING (", "").replace(")", "").split(", ")
        rearranged_pairs = []
        for pair in pairs:
            split_pair = pair.split()
            rearranged_pair = f"{split_pair[1]} {split_pair[0]}"
            rearranged_pairs.append(rearranged_pair)
        return ", ".join(rearranged_pairs)
    else:
        return line_string

"""
ΚΕΝΤΡΙΚΟ ΠΡΟΓΡΑΜΜΑ 
Φάση 1. Δημιουργία και επεξεργασία αρχείου Riga_Routes

"""

while True:
    try:
        file_name = input("Παρακαλώ εισάγετε το όνομα του αρχείου (χωρίς την κατάληξη .xlsx): ")
        gps_df = pd.read_excel(f"data/{file_name}.xlsx", dtype = object)
        break  # Εάν η ανάγνωση του αρχείου είναι επιτυχής, βγες από το βρόγχο
    except Exception as e:
        print("Παρακαλώ ελέγξτε το όνομα του αρχείου και προσπαθήστε ξανά.")

start_time = time.time()

place ="Riga, Latvia"
#  Γράφος με τα οδικά τμήματα της Ρήγας
G = ox.graph_from_place(place, network_type='drive', which_result=1)
edges_data = ox.graph_to_gdfs(G, nodes=False)

# Επιλογή συγκεκριμένων στηλών από τα δεδομένα των ακμών
columns_to_keep = ['name', "maxspeed", 'length', 'geometry']
edges = edges_data[columns_to_keep]

edges.to_excel('Riga_Routes.xlsx', index=False)
excel_file = 'Riga_Routes.xlsx'
roads_df = pd.read_excel(excel_file)

for column in roads_df.columns:
    roads_df[column] = roads_df[column].apply(rearrange_coordinates)

roads_df['pairs'] = roads_df['geometry'].apply(lambda x: [tuple(map(float, pair.split())) for pair in x.split(', ')])

# Εύρεση μέγιστων και ελαχίστων τιμών x (Γεωγ.Μήκος) και y (Γεωγ.Πλάτος) από το geometry

#GCN
roads_df['y_max'] = roads_df['pairs'].apply(lambda pairs: max(pair[0] for pair in pairs)) + 0.0018 
#GCE
roads_df['x_max'] = roads_df['pairs'].apply(lambda pairs: max(pair[1] for pair in pairs)) + 0.0034 
#GCS
roads_df['y_min'] = roads_df['pairs'].apply(lambda pairs: min(pair[0] for pair in pairs)) - 0.0018 
#GCW
roads_df['x_min'] = roads_df['pairs'].apply(lambda pairs: min(pair[1] for pair in pairs)) - 0.0034 

roads_df['geometry'] = roads_df['geometry'].str.replace(',', '')

roads_df.drop(columns=['pairs'], inplace=True)

# Προσθήκη στήλης osmID με τιμή αύξοντα αριθμό που ξεκινάει από το 1
roads_df.insert(0, 'osmID', range(1, len(roads_df) + 1))

# Αποθήκευση του DataFrame στο Excel
roads_df.to_excel('Riga_Routes.xlsx', index=False)

elapsed_time1 = time.time() - start_time
print("Η Φάση 1 εκτελέστηκε σε", elapsed_time1, "δευτερόλεπτα.")


"""
Φάση 2. Δημιουργία και επεξεργασία αρχείου GPS_Riga_Data
Διαγραφή GPS στιγμάτων 

"""

# Αρχικοποίηση στηλών και ταξινόμηση δεδομένων αρχείου GPS_Riga_Data

start_time2 = time.time()
 
gps_df['DistZ'] = None   
gps_df['Pre_CandidatesID'] = None
gps_df['CandidatesID'] = None 
gps_df['CandidatesProj'] = None
gps_df['CandidatesPGC'] = None
gps_df['CandidatesEmitProb'] = None
gps_df['CandidatesTransProb'] = None
gps_df['SelRoad'] = None
gps_df['SelNode'] = None
gps_df.at[0, 'SelNode'] = file_name
sorted_df = gps_df.sort_values(by=['VehicleID', 'SentDate', 'VehicleMessageID'], inplace=True)

prin =len(gps_df)


# Διαγραφή των GPS στιγμάτων όπου WGS84Fi ή WGS84La είναι μηδέν

gps_df = gps_df[(gps_df['WGS84Fi'] != 0) & (gps_df['WGS84La'] != 0)]
print("οι εγγραφες με WGS84Fi ή WGS84La = 0 είναι:" , int(prin) - int(len(gps_df)))

# τυπική απόκλιση
sigma_z = 4.07

# Δημιουργία R-tree index
idx = index.Index()

# Εισαγωγή των δρόμων στο R-tree index
for i, road in roads_df.iterrows():
    idx.insert(i, (road['x_min'], road['y_min'], road['x_max'], road['y_max']))


# Εύρεση των προ-υποψήφιων οδικών τμημάτων για κάθε στίγμα GPS

for index_gps, row_gps in gps_df.iterrows():
    WGS84Fi = row_gps['WGS84Fi']
    WGS84La = row_gps['WGS84La']
    candidates = list(idx.intersection((WGS84La, WGS84Fi, WGS84La, WGS84Fi)))
    for i in range(len(candidates)):
        candidates[i] += 1
    pre_candidates_ids_str = ''
    for candidate in sorted(candidates):
        pre_candidates_ids_str += str(candidate) + ' '

    pre_candidates_ids_str = pre_candidates_ids_str.strip()
    # Αντικατάσταση της τιμής της στήλης Pre_CandidatesID στο DataFrame gps_df με το pre_candidates_ids_str
    gps_df.at[index_gps, 'Pre_CandidatesID'] = pre_candidates_ids_str

gps_df.to_excel('GPS_Data.xlsx', index=False)
gps_df = pd.read_excel('GPS_Data.xlsx', dtype = object)
prin = len(gps_df)
empty_pre_candidates = []

for index, row in gps_df.iterrows():
    if pd.isna(row['Pre_CandidatesID']):
        empty_pre_candidates.append(index)


# Διαγραφή των στιγμάτων GPS χωρίς προ-υποψήφια οδικά τμήματα

gps_df.drop(empty_pre_candidates, inplace=True)
print("οι εγγραφες με κενο το Pre_CandidatesID είναι:" , int(prin) - int(len(gps_df)))


# Υπολογισμός great circle απόστασης DistZ των διαδοχικών GPS στιγμάτων |Zt-Zt+1| ανά όχημα και
# διαγραφή των GPS στιγμάτων με (DistZ >= 0 και DistZ < 2*σ_z) ή DistZ που αντιστοιχεί σε ταχύτητα οχήματος > 180km/h

row = 0
row1 = 0
count_gps = len(gps_df)
while row < count_gps :
    i = 1
    while row + i < count_gps and (gps_df.iloc[row, 2] == gps_df.iloc[row + i, 2]):
        lat1 = gps_df.iloc[row, 3]
        lon1 = gps_df.iloc[row, 4]
        lat2 = gps_df.iloc[row + i, 3]
        lon2 = gps_df.iloc[row + i, 4]
        distance = great_circle_distance(lat1, lon1, lat2, lon2)
        if distance >= 0 and distance < 2 * sigma_z:
            row1 = row + i
            gps_df.iloc[row, 6] = ''
            gps_df.iloc[row1, 6] = 1
            i = i + 1
        else:
            sent_date1 = pd.to_datetime(gps_df.iloc[row, 1], dayfirst=True)
            sent_date2 = pd.to_datetime(gps_df.iloc[row + i, 1], dayfirst=True)

            time_diff_seconds = (sent_date2 - sent_date1).total_seconds()
            max_distance = time_diff_seconds * 50
            if distance < max_distance:
                gps_df.iloc[row, 6] = distance
                row = row + i
                i = 1
            else:
                row1 = row + i
                gps_df.iloc[row, 6] = ''
                gps_df.iloc[row1, 6] = 1
                i = i + 1
    if row1 > row:
        row = row1
    row = row + 1

prin = len(gps_df)
gps_df.drop(gps_df[(gps_df['DistZ'] == 1)].index, inplace=True)
print("οι εγγραφες με (DistZ >= 0 και DistZ < 2*σ_z) ή ταχύτητα οχήματος > 180km/h  είναι:", int(prin) - int(len(gps_df)))


# Καταχώρηση στο excel αρχείο GPS_Data

gps_df.to_excel('GPS_Data.xlsx', index=False)
elapsed_time2 = time.time() - start_time2
print("Η Φάση 2 εκτελέστηκε σε", elapsed_time2, "δευτερόλεπτα.")
print("Το συνολικό πρόγραμμα εκτελέστηκε σε", elapsed_time1 + elapsed_time2, "δευτερόλεπτα.")

