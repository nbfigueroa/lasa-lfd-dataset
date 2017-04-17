# lasa-lfd-dataset
Kinesthetic Recordings of Complex Sequential Manipulation Tasks. This dataset is intended to test and benchmark segmentation, clustering and learning algorithms.

---

## Complex Sequential Tasks
The dataset currently contains time-series of kinesthetic demonstrations performed at the EPFL-LASA lab with KUKA-LWR robots and f/t sensors. 

### Tasks
We currently have recordings of 3 complex sequential tasks that have been segmented, processed and learned.
- Vegetable Grating
- Dough Rolling
- Zucchini Peeling

Each task has it's own folder, whose structure is as follows:
```
.
└── continous
└── per-task
└── per-action
```

- ```continous```: Here you will find an unsegmented Continous Complex Sequential Task Recording of a complete high-level task.
From initial to final goal/configuration. For example, multiple iterations of rolling sequences towards the desired dough diamter or multiple peeling sequences until the whole zucchini is peeled.
- ```per-task```: Here I provide pre-segmented sequences tasks, i.e. each time-series is a sequence of reach-grate-back for the grating demonstrations.
- ```per-action```: Here you will find fully segmented time-series corresponding to the same action within a task

### Observed Data

#### Vegetable Grating 
Demonstration of a Carrot Grating Task consisting of 12 (7-d) time-series X = {x_1,..,x_T} with variable length T.  
**Dimensions:** x = {pos_x, pos_y, pos_z, q_i, q_j, q_k, q_w}


#### Dough Rolling
Demonstration of a Dough Rolling Task consisting of 15 (12-d) time-series X = {x_1,..,x_T} with variable length T.  
**Dimensions:** x = {pos_x, pos_y, pos_z, roll, pitch , yaw, f_x, f_y, f_z, tau_x, tau_y, tau_z}

#### Zucchini Peeling
Demonstration of a Peeling Task consisting of 2 (13-d) time-series X = {x_1,..,x_T} with variable length T.  
**Dimensions:** x = {pos_x, pos_y, pos_z, q_i, q_j, q_k, q_w, f_x, f_y, f_z, tau_x, tau_y, tau_z}

---

## Task-Wrench-Space Primitives
The dataset currently contains Task-Wrench-Space Primitives, which are 6-dimensional Ellipsoids indicating force/torque space used in specific tasks.

### Tasks
- Cuting
- Screwing
- Drawing
