# Robot Manipulators for Sterilisation Control Tests
A 6-DOF robotic prototype automating the sterilisation control test with biological indicators, aimed at small- and medium-sized hospitals.

## Demo


https://github.com/user-attachments/assets/21a6fb4d-3214-4c9d-89f8-d67e49fbd7de


## Problem
Sterilisation technicians handle large volumes of repetitive pick-and-place tasks. These include sorting instruments, moving boxes, running control tests etc.. Studies show that 88.3% of sterilisation department staff experience musculoskeletal disorders as a result. Existing automated solutions (e.g. Gibotech's system) solve this, but their cost, size, and complexity make them inaccessible to small- and medium-sized hospitals.

**Problem formulation:** How can a robot manipulator be used to develop a prototype for performing a sterilisation control test with biological indicators in small- to medium-sized hospitals?

## Key Features
- Automated insertion of biological indicator (BI) samples into incubator slots at a fixed 45° angle.
- Automated placement and removal of the incubator lid.
- Automated retrieval of incubated samples into an inspection tray.
- Custom v-shaped parallel gripper (135°) for self-centring cylindrical samples.
- Trajectory planning using Linear Segments with Parabolic Blends (LSPB).
- Forward and inverse kinematics model, verified against the robot's encoder data.

## Tech Stack
- **AgileX Piper** — 6-DOF robotic manipulator.
- **Python (PiperSDK)** — robot control interface.
- **Python (Robotics Toolbox)** - Simulating and valdiating Forwars & Inverse kinematics.
- **Python (NumPy) & (time)** - Plan, Caluclate and send positions to the robot for trajectory planning.
- **CAN bus** — robot-to-PC communication.
- **Custom 3D-printed parallel gripper** - v-shaped jaws for self-centering of BI's.
  
## My Role & Skills Demonstrated
- Deriving Forward and Inverse kinematics and verifiying them.
- Creating Trajectory planning code based on given formulas.
- Conducting testing of trajectory planning.
- Creating the gripper control code. 

## How It Works
Joint positions for each target pose were calculated offline using inverse kinematics and hardcoded into the control script, together with a count (movement duration) and speed for each motion. For multi-point movements, an LSPB (Linear Segments with Parabolic Blends) function takes a start pose, end pose, segment duration, and maximum acceleration, and generates the trajectory sent to the robot. The robot's own onboard forward kinematics are used to confirm it reached the intended pose. The gripper is triggered manually (button press) in sync with each automated movement.

Sequence performed by the prototype:
1. Pick a BI sample from the sample-holding block
2. Crush BI sample and insert into incubator slot
3. Place the incubator lid onto the incubator
4. Wait through the incubation cycle (~48 hours)
5. Remove the lid and return it to the lid holder
6. Retrieve the incubated samples and place them in the inspection tray

## Setup

**Requirements**
- Ubuntu (native install — dual boot or bare metal; the CAN interface needs direct hardware access)
- AgileX Piper 6-DOF manipulator + CAN-to-USB adapter
- Python 3 + pip
- Custom 3D-printed parallel gripper

**Software installation**
```bash
# Install the Piper SDK
pip3 install piper_sdk

# Install CAN tools
sudo apt update && sudo apt install can-utils ethtool

# Activate the CAN interface (bitrate 1,000,000 is the Piper default)
bash can_activate.sh can0 1000000
```

**Physical setup**
1. Mount the robot base to a stable, fixed surface.
2. Mount the gripper onto the robot's wrist mounting point.
3. Connect the CAN-to-USB adapter (CAN_H, CAN_L, VCC, GND) between the robot and the PC.
4. Position the workspace fixtures: BI-holding block, incubator, incubator lid holder, and BI inspection tray.

## Usage
[ASK: what command(s) does someone run to execute a cycle?]

## Results
| Test | Outcome |
|---|---|
| Rotation Test | Passed — ±0.24 mm (x-axis) / ±0.40 mm (z-axis) deviation |
| Linear Path Test | Failed — up to 5.6 mm deviation (joint-space plan ≠ straight line in Cartesian space) |
| Payload Test | Passed |
| Stability Test | Passed — BI and lid, 20–100% speed |
| Complete Solution Test | Passed — 39/40 trials (1 failure: gripper missed sample retrieval) |


## Challenges & What I Learned
- Linear vertical trajectories planned in joint space did not translate to straight-line motion in Cartesian space — the rotation test (planned in Cartesian space) performed better, confirming this.
- Long, slim gripper jaws (needed to avoid collisions between samples) caused shaking and bending under load.
- The incubator's height limited the manipulator's reachable orientations in parts of the workspace.
- One of the team members didn't put enough effort into the project in the first half of the semester. This was solved by confronting the problem and laying out a plan for said team member. This improved our efficiency as a group.
- Workspece constraints presented difficulties for the execution of some of the insertations of BI samples into the incubator slots. Simulating the workspace with the robots joint constraints before implementing the solution would resolve this issue.
  

## Future Improvements
- Automate the full sterilisation control test, which includes automating the gripper sequence and workspace issue. 
- Redesign gripper jaws with a more durable material to reduce bending. 
- Integrate gripper control into the main control software (currently a separate manual button).
- Add an emergency stop button.
- Switch trajectory planning to Cartesian space for improved accuracy.
- Simulate the workspace and configure the positions of BI-holder, Incubaotr and BI-tray for better results.
- Include autormation of retrieval of BI samples directly from the autoclave machine. 

## License & Contact
[ASK: license? list all 5 teammates + emails, or just you?]
