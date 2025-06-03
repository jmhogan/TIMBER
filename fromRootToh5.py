import ROOT
from ROOT import TFile
import uproot
import sys, os
import gc
import h5py
import awkward as ak
import numpy as np
from sklearn.model_selection import train_test_split

def awkward_to_fixed_numpy(arr, length=16, fill_value=-1):
    """
    Convert a 2D awkward array to a NumPy array with fixed-length rows.

    Parameters:
        arr (ak.Array): A 2D awkward array with variable-length inner lists.
        length (int): The target fixed length of each row.
        fill_value (int or float): Value to use for padding.

    Returns:
        np.ndarray: A (N, length) NumPy array.
    """
    # Pad to fixed length (with None)
    padded = ak.pad_none(arr, length, clip=True)
    # Replace None with fill_value
    filled = ak.fill_none(padded, fill_value)
    # Convert to NumPy
    return ak.to_numpy(filled)


# File to convert 
# RDF_BprimeBprimeto2B4Tau_MB-1000_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2022_0.root 
# into 3 .h5 files for Topograph test

# List of filenames
filenames = [
    "RDF_BprimeBprimeto2B4Tau_MB-1000_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2022_0.root",
    "RDF_BprimeBprimeto2B4Tau_MB-1000_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2022EE_0.root",
    "RDF_BprimeBprimeto2B4Tau_MB-1000_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023_0.root",
    "RDF_BprimeBprimeto2B4Tau_MB-1600_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023_41.root",
    "RDF_BprimeBprimeto2B4Tau_MB-400_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023_61.root",
    #"RDF_BprimeBprimeto2B4Tau_MB-1600_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023_71.root", This one errored out :(
    "RDF_BprimeBprimeto2B4Tau_MB-700_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023_75.root",
    "RDF_BprimeBprimeto2B4Tau_MB-1000_MXi-2000_TuneCP5_13p6TeV-madgraph-pythia8_2023BPix_0.root",
]

# Initialize lists
ObjectList_indicesAK = []
matchabilityArr = []
nObjects = []
nbjets = []
matchability = []
PtListObject = []
EtaListObject = []
PhiListObject = []
EnergyListObject = []
TaggedListObject = []
pdgIdListGen = []
ptListGen = []
etaListGen = []
phiListGen = []
massListGen = []

for filename in filenames:
    print(f"Opening...{filename}")
    f = uproot.open(filename)
    events = f['Events_Nominal;1']
    nevents = events.num_entries
    print(f"{nevents = }")

    # Extend master lists with data from current file
    ObjectList_indicesAK.extend(events['ObjectList_indices'].array())
    nObjects.extend(events['nObjects'].array())
    nbjets.extend(events['nbjets'].array())
    matchability.extend(events['matchability'].array())
    PtListObject.extend(events['PtListObject'].array())
    EtaListObject.extend(events['EtaListObject'].array())
    PhiListObject.extend(events['PhiListObject'].array())
    EnergyListObject.extend(events['EnergyListObject'].array())
    TaggedListObject.extend(events['TaggedListObject'].array())
    pdgIdListGen.extend(events['pdgIdListGen'].array())
    ptListGen.extend(events['ptListGen'].array())
    etaListGen.extend(events['etaListGen'].array())
    phiListGen.extend(events['phiListGen'].array())
    massListGen.extend(events['massListGen'].array())

# Now each variable contains the concatenated arrays from all files!
print(f"Total events: {len(nObjects)}")
numOfEvents = len(nObjects)


ObjectList_indices = awkward_to_fixed_numpy(ObjectList_indicesAK)

Objects = np.zeros(shape=(len(PtListObject), 16, 5))
for i in range(0, len(PtListObject)): # 27 at the moment
    event = np.zeros(shape=(16, 5))
    for j in range(0, len(PtListObject[i])):
        newObject = np.array([PtListObject[i][j], EtaListObject[i][j], PhiListObject[i][j], EnergyListObject[i][j], TaggedListObject[i][j]])
        event[j] = newObject
    for k in range(len(PtListObject[i]), 16):
        # Pad with 0s (Checked the other Topograph input)
        newObject = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        event[k] = newObject
    Objects[i] = event
    
    
partons = np.zeros(shape=(len(pdgIdListGen), 10, 5))
for i in range(0, len(pdgIdListGen)): # 27 at the moment
    event = np.zeros(shape=(10, 5))
    for j in range(0, len(pdgIdListGen[i])):
        newParton = np.array([pdgIdListGen[i][j], ptListGen[i][j], etaListGen[i][j], phiListGen[i][j], massListGen[i][j]])
        event[j] = newParton
    for k in range(len(pdgIdListGen[i]), 10):
        newParton = np.array([0.0, 0.0, 0.0, 0.0, 0.0])
        event[k] = newParton
    partons[i] = event


# Split events before passing into the three h5 files

from sklearn.model_selection import train_test_split

# Step 1: Convert nObjects, matchability, and nbjets to NumPy arrays and flatten
nObjects_np = np.array(nObjects).flatten()
matchability_np = np.array(matchability).flatten()
nbjets_np = np.array(nbjets).flatten()

# Step 2: Filter rare matchability values
value_counts = np.unique(matchability_np, return_counts=True)
valid_values = value_counts[0][value_counts[1] >= 5]
valid_mask = np.isin(matchability_np, valid_values)
nObjects_np = nObjects_np[valid_mask]
matchability_np = matchability_np[valid_mask]
nbjets_np = nbjets_np[valid_mask]
ObjectList_indices = ObjectList_indices[valid_mask]
Objects = Objects[valid_mask]
partons = partons[valid_mask]
numOfEvents = len(nObjects_np)
print(f"Filtered {np.sum(~valid_mask)} events with rare matchability values (<5 occurrences)")
print("Filtered matchability distribution:", np.unique(matchability_np, return_counts=True))

# Step 3: Create stratification keys
# Bin matchability for stratification (finer bins)
matchability_binned = np.digitize(matchability_np, bins=np.linspace(0, 100, 7))  # 6 bins
# Combine for stratification
stratify_key = np.vstack([matchability_binned]).T

# Step 4: Perform stratified splitting
try:
    # Split into train (60%) and val+test (40%)
    train_idx, val_test_idx = train_test_split(
        np.arange(numOfEvents),
        test_size=0.4,
        stratify=stratify_key,
        random_state=42
    )
    # Split val+test into val (20%) and test (20%)
    val_idx, test_idx = train_test_split(
        val_test_idx,
        test_size=0.5,
        stratify=stratify_key[val_test_idx],
        random_state=42
    )
except ValueError as e:
    print(f"Stratified splitting failed: {e}. Falling back to random splitting.")
    train_idx, val_test_idx = train_test_split(
        np.arange(numOfEvents),
        test_size=0.4,
        random_state=42
    )
    val_idx, test_idx = train_test_split(
        val_test_idx,
        test_size=0.5,
        random_state=42
    )

# Step 5: Balance matchability value 63
value_63_mask = matchability_np == 63
value_63_indices = np.where(value_63_mask)[0]
other_indices = np.where(~value_63_mask)[0]
train_63_indices, val_test_63_indices = train_test_split(
    value_63_indices, test_size=0.4, random_state=42
)
val_63_indices, test_63_indices = train_test_split(
    val_test_63_indices, test_size=0.5, random_state=42
)
train_other_indices, val_test_other_indices = train_test_split(
    other_indices, test_size=0.4, random_state=42
)
val_other_indices, test_other_indices = train_test_split(
    val_test_other_indices, test_size=0.5, random_state=42
)
# Combine indices
train_idx = np.concatenate([train_63_indices, train_other_indices])
val_idx = np.concatenate([val_63_indices, val_other_indices])
test_idx = np.concatenate([test_63_indices, test_other_indices])

# Step 6: Initialize arrays for each split
# Training set
ObjectList_indicesTrain = ObjectList_indices[train_idx]
matchabilityTrain = matchability_np[train_idx]
nObjectsTrain = nObjects_np[train_idx]
nbjetsTrain = nbjets_np[train_idx]
ObjectsTrain = Objects[train_idx]
partonsTrain = partons[train_idx]

# Validation set
ObjectList_indicesVal = ObjectList_indices[val_idx]
matchabilityVal = matchability_np[val_idx]
nObjectsVal = nObjects_np[val_idx]
nbjetsVal = nbjets_np[val_idx]
ObjectsVal = Objects[val_idx]
partonsVal = partons[val_idx]

# Test set
ObjectList_indicesTest = ObjectList_indices[test_idx]
matchabilityTest = matchability_np[test_idx]
nObjectsTest = nObjects_np[test_idx]
nbjetsTest = nbjets_np[test_idx]
ObjectsTest = Objects[test_idx]
partonsTest = partons[test_idx]

# Step 7: Log distributions to verify balance
print(f"Training set size: {len(train_idx)} ({len(train_idx)/numOfEvents:.2%})")
print(f"Validation set size: {len(val_idx)} ({len(val_idx)/numOfEvents:.2%})")
print(f"Test set size: {len(test_idx)} ({len(test_idx)/numOfEvents:.2%})")
print("\nDistribution of njets:")
print(f"Training njets mean: {np.mean(nObjectsTrain):.2f}, std: {np.std(nObjectsTrain):.2f}")
print(f"Validation njets mean: {np.mean(nObjectsVal):.2f}, std: {np.std(nObjectsVal):.2f}")
print(f"Test njets mean: {np.mean(nObjectsTest):.2f}, std: {np.std(nObjectsTest):.2f}")
print("\nDistribution of matchability:")
train_match_counts = np.unique(matchabilityTrain, return_counts=True)
val_match_counts = np.unique(matchabilityVal, return_counts=True)
test_match_counts = np.unique(matchabilityTest, return_counts=True)
print(f"Training matchability: {train_match_counts}")
print(f"Validation matchability: {val_match_counts}")
print(f"Test matchability: {test_match_counts}")
print(f"Training matchability 63 proportion: {np.sum(matchabilityTrain == 63) / len(matchabilityTrain):.4f}")
print(f"Validation matchability 63 proportion: {np.sum(matchabilityVal == 63) / len(matchabilityVal):.4f}")
print(f"Test matchability 63 proportion: {np.sum(matchabilityTest == 63) / len(matchabilityTest):.4f}")
print("\nDistribution of nbjets:")
print(f"Training nbjets mean: {np.mean(nbjetsTrain):.2f}, std: {np.std(nbjetsTrain):.2f}")
print(f"Validation nbjets mean: {np.mean(nbjetsVal):.2f}, std: {np.std(nbjetsVal):.2f}")
print(f"Test nbjets mean: {np.mean(nbjetsTest):.2f}, std: {np.std(nbjetsTest):.2f}")

# --- Define structured dtypes ---
jet_dtype = np.dtype([
    ('pt', 'f4'),
    ('eta', 'f4'),
    ('phi', 'f4'),
    ('energy', 'f4'),
    ('is_tagged', '?')
])

parton_dtype = np.dtype([
    ('PDGID', 'i4'),
    ('pt', 'f4'),
    ('eta', 'f4'),
    ('phi', 'f4'),
    ('mass', 'f4')
])

# --- Conversion functions ---
def convert_to_structured_jets(jet_array):
    if jet_array.shape[0] == 0:
        return np.zeros((0, jet_array.shape[1]), dtype=jet_dtype)
    structured = np.zeros(jet_array.shape[:2], dtype=jet_dtype)
    structured['pt'] = jet_array[:, :, 0].astype(np.float32)
    structured['eta'] = jet_array[:, :, 1].astype(np.float32)
    structured['phi'] = jet_array[:, :, 2].astype(np.float32)
    structured['energy'] = jet_array[:, :, 3].astype(np.float32)
    structured['is_tagged'] = jet_array[:, :, 4].astype(bool)
    return structured

def convert_to_structured_partons(parton_array):
    if parton_array.shape[0] == 0:
        return np.zeros((0, parton_array.shape[1]), dtype=parton_dtype)
    structured = np.zeros(parton_array.shape[:2], dtype=parton_dtype)
    structured['PDGID'] = parton_array[:, :, 0].astype(np.int32)
    structured['pt']    = parton_array[:, :, 1].astype(np.float32)
    structured['eta']   = parton_array[:, :, 2].astype(np.float32)
    structured['phi']   = parton_array[:, :, 3].astype(np.float32)
    structured['mass']  = parton_array[:, :, 4].astype(np.float32)
    return structured

# --- Convert arrays as needed ---
ObjectList_indicesTrain = ObjectList_indicesTrain.astype(np.int16)
ObjectList_indicesVal = ObjectList_indicesVal.astype(np.int16)
ObjectList_indicesTest = ObjectList_indicesTest.astype(np.int16)
matchabilityTrain = matchabilityTrain.astype(np.int16)
matchabilityVal = matchabilityVal.astype(np.int16)
matchabilityTest = matchabilityTest.astype(np.int16)
nbjetsTrain = nbjetsTrain.astype(np.int16)
nbjetsVal = nbjetsVal.astype(np.int16)
nbjetsTest = nbjetsTest.astype(np.int16)
nObjectsTrain = nObjectsTrain.astype(np.int16)
nObjectsVal = nObjectsVal.astype(np.int16)
nObjectsTest = nObjectsTest.astype(np.int16)
partonsTrain = convert_to_structured_partons(partonsTrain)
partonsVal = convert_to_structured_partons(partonsVal)
partonsTest = convert_to_structured_partons(partonsTest)
ObjectsTrain = convert_to_structured_jets(ObjectsTrain)
ObjectsVal = convert_to_structured_jets(ObjectsVal)
ObjectsTest = convert_to_structured_jets(ObjectsTest)

with h5py.File("b_train.h5", 'w') as f:
    grp = f.create_group('delphes')
    grp.create_dataset('jets', data=ObjectsTrain)
    grp.create_dataset('jets_indices', data=ObjectList_indicesTrain)
    grp.create_dataset('matchability', data=matchabilityTrain)
    grp.create_dataset('njets', data=nObjectsTrain)
    grp.create_dataset('nbjets', data=nbjetsTrain)
    grp.create_dataset('partons', data=partonsTrain)
    print("b_train.h5 file created successfully.")
    
with h5py.File("b_val.h5", 'w') as f:
    grp = f.create_group('delphes')
    grp.create_dataset('jets', data=ObjectsVal)
    grp.create_dataset('jets_indices', data=ObjectList_indicesVal)
    grp.create_dataset('matchability', data=matchabilityVal)
    grp.create_dataset('njets', data=nObjectsVal)
    grp.create_dataset('nbjets', data=nbjetsVal)
    grp.create_dataset('partons', data=partonsVal)
    print("b_val.h5 file created successfully.")
    
with h5py.File("b_test.h5", 'w') as f:
    grp = f.create_group('delphes')
    grp.create_dataset('jets', data=ObjectsTest)
    grp.create_dataset('jets_indices', data=ObjectList_indicesTest)
    grp.create_dataset('matchability', data=matchabilityTest)
    grp.create_dataset('njets', data=nObjectsTest)
    grp.create_dataset('nbjets', data=nbjetsTest)
    grp.create_dataset('partons', data=partonsTest)
    print("b_test.h5 file created successfully.")