import time
import csv
import numpy as np
from piper_sdk import *


log_counter = 1
speed = 15                #Speed parameter, vi kan bruge til at styre piperens hastighed med 
count_constant = 100

# LSPB parametre
TRAJ_STREAM_DT = 0.02        
TRAJ_STREAM_SPEED = 90      


# ---------------------------------------------------------------------- 
# LSPB Funktioner 
# ---------------------------------------------------------------------


def lspb_via_points(via, segment_durations, max_accel, dt=0.01):
    via = np.asarray(via, dtype=float)
    td = np.asarray(segment_durations, dtype=float)
    n = len(via)
    assert len(td) == n - 1, "Der skal være n-1 segment durations for n via points"

    v = np.zeros(n - 1)
    t_b = np.zeros(n)
    a_signed = np.zeros(n)

    #Første via point (start from rest)
    if abs(via[1] - via[0]) < 1e-12:
        t_b[0] = 0.0; a_signed[0] = 0.0; v[0] = 0.0
    else:
        a_signed[0] = np.sign(via[1] - via[0]) * max_accel
        disc = td[0]**2 - 2 * (via[1] - via[0]) / a_signed[0]
        if disc < 0:
            raise ValueError("Det er ikke muligt. Der er enten for lav acceleration eller for kort segmentvarighed(Første via)")
        t_b[0] = td[0] - np.sqrt(disc)
        v[0] = (via[1] - via[0]) / (td[0] - 0.5 * t_b[0])

    # sidste via point
    if abs(via[-1] - via[-2]) < 1e-12:
        t_b[-1] = 0.0; a_signed[-1] = 0.0; v[-1] = 0.0
    else:
        a_signed[-1] = np.sign(via[-2] - via[-1]) * max_accel
        disc = td[-1]**2 + 2 * (via[-1] - via[-2]) / a_signed[-1]
        if disc < 0:
            raise ValueError("Det er ikke muligt. Der er enten for lav acceleration eller for kort segmentvarighed(Første via)")
        t_b[-1] = td[-1] - np.sqrt(disc)
        v[-1] = (via[-1] - via[-2]) / (td[-1] - 0.5 * t_b[-1])

    # mellem via points
    for k in range(1, n - 1):
        v_in = v[0] if k == 1 else (via[k] - via[k - 1]) / td[k - 1]
        v_out = v[-1] if k == n - 2 else (via[k + 1] - via[k]) / td[k]
        v[k - 1] = v_in
        v[k] = v_out
        if abs(v_out - v_in) < 1e-12:
            t_b[k] = 0.0; a_signed[k] = 0.0
        else:
            a_signed[k] = np.sign(v_out - v_in) * max_accel
            t_b[k] = (v_out - v_in) / a_signed[k]

    # Linear segment varighed
    t_lin = np.zeros(n - 1)
    t_lin[0] = td[0] - t_b[0] - 0.5 * t_b[1]
    t_lin[-1] = td[-1] - t_b[-1] - 0.5 * t_b[-2]
    for j in range(1, n - 2):
        t_lin[j] = td[j] - 0.5 * t_b[j] - 0.5 * t_b[j + 1]
    if any(t_lin < -1e-9):
        raise ValueError("Det er ikke muligt. Der er enten for lav acceleration eller for kort segmentvarighed(lineær segment)")

   
    T_total = float(np.sum(td))
    t = np.arange(0, T_total + dt, dt)
    theta = np.zeros_like(t)
    theta_dot = np.zeros_like(t)
    theta_ddot = np.zeros_like(t)
    via_times = np.concatenate(([0.0], np.cumsum(td)))
    blend_start = via_times - t_b / 2
    blend_end = via_times + t_b / 2
    blend_start[0] = 0.0
    blend_end[-1] = T_total

    for i, ti in enumerate(t):
        in_blend = False
        for k in range(n):
            if blend_start[k] <= ti <= blend_end[k]:
                in_blend = True
                tau = ti - blend_start[k]
                v_in = 0.0 if k == 0 else v[k - 1]
                if k == 0:
                    theta_at_start = via[0]
                else:
                    theta_at_start = (via[k] - v_in * t_b[k] / 2
                                      - 0.5 * a_signed[k] * (t_b[k] / 2)**2)
                theta[i] = theta_at_start + v_in * tau + 0.5 * a_signed[k] * tau**2
                theta_dot[i] = v_in + a_signed[k] * tau
                theta_ddot[i] = a_signed[k]
                break
        if not in_blend:
            for k in range(n - 1):
                if blend_end[k] < ti < blend_start[k + 1]:
                    tau = ti - blend_end[k]
                    if k == 0:
                        theta_at_blend_end = via[0] + 0.5 * a_signed[0] * t_b[0]**2
                    else:
                        v_in = v[k - 1]
                        theta_at_blend_end = (via[k] + v_in * t_b[k] / 2
                                              + 0.5 * a_signed[k] * (t_b[k] / 2)**2)
                    theta[i] = theta_at_blend_end + v[k] * tau
                    theta_dot[i] = v[k]
                    theta_ddot[i] = 0.0
                    break

    return t, theta, theta_dot, theta_ddot, t_b


    # ----------------------------------------------------------------------------
    # MAIN PROGRAM 
    # ---------------------------------------------------------------------------

if __name__ == "__main__":
    piper = C_PiperInterface_V2("can0")
    piper.ConnectPort()
    while not piper.EnablePiper():
        time.sleep(0.01)

    factorRadToDeg = 57295.7795       # 1000*180/pi  (rad → milli-degrees)
    factorMeterToMicro = 1000000      # m → μm
    position = [0, 0, 0, 0, 0, 0]
    count = 0

    # ------------------------------------------------------------------------
    # CSV setup
    # ------------------------------------------------------------------------
    csv_filename = "Fuld_Loesning_Koordinater"
    with open(csv_filename, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow([
            "ID", "q1_rad", "q2_rad", "q3_rad", "q4_rad", "q5_rad", "q6_rad",
            "x", "y", "z", "Rx", "Ry", "Rz"
        ])
    print("CSV-fil oprettet:", csv_filename)

    def log_to_csv(q, pose):
        global log_counter
        with open(csv_filename, mode="a", newline="") as file:
            writer = csv.writer(file)
            writer.writerow([
                log_counter,
                round(q[0], 6), round(q[1], 6), round(q[2], 6),
                round(q[3], 6), round(q[4], 6), round(q[5], 6),
                round(pose[0], 6), round(pose[1], 6), round(pose[2], 6),
                round(pose[3], 6), round(pose[4], 6), round(pose[5], 6)
            ])
        log_counter += 1

    # ------------------------------------------------------------------------
    # Joint movement code fra PiperSDK GIthub
    # ------------------------------------------------------------------------
    def piper_joint_movement(jointRadians):
        joint_0 = round(jointRadians[0] * factorRadToDeg)
        joint_1 = round(jointRadians[1] * factorRadToDeg)
        joint_2 = round(jointRadians[2] * factorRadToDeg)
        joint_3 = round(jointRadians[3] * factorRadToDeg)
        joint_4 = round(jointRadians[4] * factorRadToDeg)
        joint_5 = round(jointRadians[5] * factorRadToDeg)
        print(joint_0, joint_1, joint_2, joint_3, joint_4, joint_5)
        piper.MotionCtrl_2(0x01, 0x01, speed, 0x00)
        piper.JointCtrl(joint_0, joint_1, joint_2, joint_3, joint_4, joint_5)



    # ------------------------------------------------------------------------
    # LSPB streaming funktioner(Disse bruger vi til at sende instrukserne til robotten)
    # ------------------------------------------------------------------------
    def execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                dt=TRAJ_STREAM_DT):
        """
        Denne del af koden laver trajectory gennem via_points og sender det til Piperen på baggrund af noget input:

        via_points:        Dette er de forskellige viapoints, der skal vælges i et array
        segment_durations: varighederne af hvert egment i et array
        max_accel:         max-accelerationen for robotten. DEtte har vi valgt ud fra det langsomste joint
        dt:                Dette er, hvor often koden sendes til robotten: 0,02s = 50 Hz 
        """
        via = np.asarray(via_points, dtype=float)
        n_joints = via.shape[1]
        assert n_joints == 6, "Der forvnetes 6-joint trajectory"

        # Normaliser den valgte max_accel til alle joints
        if np.isscalar(max_accel):
            max_accels = [float(max_accel)] * n_joints
        else:
            max_accels = list(max_accel)

        # Genrering af position trajectory for hver joint
        joint_trajs = []
        t_array = None
        for j in range(n_joints):
            t_arr, theta, _, _, _ = lspb_via_points(
                via[:, j], segment_durations, max_accels[j], dt=dt
            )
            joint_trajs.append(theta)
            if t_array is None:
                t_array = t_arr

        n_samples = len(t_array)
        print(f"  → Streaming {n_samples} samples over {t_array[-1]:.2f} s "
              f"({1/dt:.0f} Hz)")

        piper.MotionCtrl_2(0x01, 0x01, TRAJ_STREAM_SPEED, 0x00)

        # Kør igennem de forskellige points.
        loop_start = time.time()
        for i in range(n_samples):
            j0 = round(joint_trajs[0][i] * factorRadToDeg)
            j1 = round(joint_trajs[1][i] * factorRadToDeg)
            j2 = round(joint_trajs[2][i] * factorRadToDeg)
            j3 = round(joint_trajs[3][i] * factorRadToDeg)
            j4 = round(joint_trajs[4][i] * factorRadToDeg)
            j5 = round(joint_trajs[5][i] * factorRadToDeg)
            piper.JointCtrl(j0, j1, j2, j3, j4, j5)

            # Sleep indtil næste sample
            next_target = loop_start + t_array[i] + dt
            sleep_time = next_target - time.time()
            if sleep_time > 0:
                time.sleep(sleep_time)

        elapsed = time.time() - loop_start
        print(f"  → Trajectory complete in {elapsed:.2f} s")

    # ------------------------------------------------------------------------
    #     MAIN LOOP 
    # ------------------------------------------------------------------------
    while True:
        count += 1
        log_now = False

        if count == count_constant * 1:
            print("Pos 1 -----------")
            position = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #0 position (Start)
            piper_joint_movement(position)
            log_now = True

        elif count == count_constant * 2:
            print("Pos 2 -----------")
            position = [0.0, 0.9139, -0.6589, 0.0, -0.1678, 0.0]   # Home Position
            piper_joint_movement(position)
            log_now = True

        elif count == count_constant * 3:
            print("Pos 3 -----------")
            position = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            piper_joint_movement(position)
            log_now = True

        # --------------------------------------------------------------------
        #  TRAJ 1  —  Home BI holder → Grip A3 BI Holder via LSPB
        elif count == count_constant * 4:
            print("TRAJ 1 -----------  (LSPB)")

            pos_start = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            pos_end   = [0.7134, 1.6627, -0.9093, 0.0, 0.9047, 2.356] #A3 Grip position

            
            via_points = np.array([
                pos_start,
                pos_end,
            ])

            segment_durations = [3.0]   # tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
         # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 2  —  Grip A3 BI holder → Lift A3 

        elif count == count_constant * 5:
            print("TRAJ 2 -----------  (LSPB )")

            pos_start = [0.7134, 1.6627, -0.9093, 0.0, 0.9047, 2.356] #A3 Grip position
            pos_end   = [0.7069, 1.6025, -0.9841, 0.0, 1.0397, 2.356] #A3 Approach 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

         
            segment_durations = [3.0]   # tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        elif count == count_constant * 6:
            print("Pos 6 -----------")
            position = [0.1404, 1.2, -0.82, 0.0, 1.305, 1.593] #B10 approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 3  —  B10 Approach → B10 Insert sample

        elif count == count_constant * 7:  
            print("TRAJ 3 -----------  (LSPB )")

            
            pos_start = [0.1404, 1.2, -0.82, 0.0, 1.305, 2.356] #B10 Approach
            pos_end   = [0.1395, 1.2639, -0.7622, 0.0, 1.1563, 2.356] #B10 Insert sample 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            
            segment_durations = [3.0]    #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

         # --------------------------------------------------------------------
        #  TRAJ 4  —  B10 Insert sample → B10 Approach

        elif count == count_constant * 8: 
            print("TRAJ 4 -----------  (LSPB )")

            pos_start = [0.1404, 1.2621, -0.7607, 0.0, 1.1567, 2.356] #B10 Insert sample
            pos_end   = [0.1049, 1.2343, -0.8847, 0.0, 1.3084, 2.356] #B10  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            
            segment_durations = [3.0]    # tid af bevægelsen fra A til B

            max_accel = 3.14            # rad/s²

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True


        # --------------------------------------------------------------------

        elif count == count_constant * 9:
            print("Pos 3 -----------")
            position = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            piper_joint_movement(position)
            log_now = True

        # --------------------------------------------------------------------
        #  TRAJ 5  —  Home BI holder → Grip A2 BI Holder via LSPB 
        elif count == count_constant * 10: 
            print("TRAJ 1 -----------  (LSPB )")

           
            pos_start = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            pos_end   = [0.773, 1.7114, -0.9743, 0.0, 0.921, 2.402] #A2 Grip position

          
            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   # tid af bevægelsen fra A til B

            max_accel = 3.14         

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
         # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 6  —  Grip A2 BI holder → Lift A2 

        elif count == count_constant * 11:
            print("TRAJ 2 -----------  (LSPB )")

      
            pos_start = [0.7718, 1.7116, -0.9746, 0.0, 0.9211, 2.402] #A2 Grip position
            pos_end   = [0.7718, 1.6547, -1.0489, 0.0, 1.0523, 2.402] #A2 Approach 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

         
            segment_durations = [3.0]   # tid af bevægelsen fra A til B

            max_accel = 3.14          

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 12:
            print("Pos 6 -----------")
            position = [-0.1377, 1.2894, -0.8883, 0.0, 1.2569, -0.1526] #B5 approach
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 3  —  B5 Approach → B5 Insert sample

        elif count == count_constant * 13:
            print("TRAJ 3 -----------  (LSPB )")

          
            pos_start = [-0.1377, 1.2894, -0.8883, 0.0, 1.2569, -0.1526] #B5 Approach
            pos_end   = [-0.1526, 1.3404, -0.8265, 0.0, 1.1441, -0.1526] #B5 Insert sample 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

         # --------------------------------------------------------------------
        #  TRAJ 7  —  B5 Insert sample → B5 Approach

        elif count == count_constant * 14:
            print("TRAJ 4 -----------  (LSPB)")

         
            pos_start = [-0.1526, 1.3404, -0.8265, 0.0, 1.1441, -0.1526] #B5 Insert sample
            pos_end   = [-0.1526, 1.3038, -0.9511, 0.0, 1.3053, -0.1526] #B5  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

          
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14          

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 15:
            print("Pos 3 -----------")
            position = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            piper_joint_movement(position)
            log_now = True


        # --------------------------------------------------------------------

        #  TRAJ 8  —  Home BI holder → Grip A1 BI Holder via LSPB 
        elif count == count_constant * 16:
            print("TRAJ 1 -----------  (LSPB )")

            pos_start = [0.8099, 1.6641, -1.1255, 0.0, 1.1195, 0.8099]  # Home BI holder
            pos_end   = [0.8339, 1.7684, -1.0541, 0.0, 0.9437, 2.402] #A1 Grip position

       
            via_points = np.array([
                pos_start,
                pos_end,
            ])


            segment_durations = [3.0]   # tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
         # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 9  —  Grip A1 BI holder → Lift A1 

        elif count == count_constant * 17:
            print("TRAJ 2 -----------  (LSPB streaming)")

          
            pos_start = [0.8312, 1.7666, -1.0514, 0.0, 0.9429, 2.402] #A1 Grip position
            pos_end   = [0.8312, 1.7126, -1.1245, 0.0, 1.07, 2.402] #A1 Approach 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        elif count == count_constant * 18:
            print("Pos 6 -----------")
            position = [0.0222, 1.4202, -1.0577, 0.0, 1.2955, 1.593] #B1 approach
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 10  —  B1 Approach → B1 Insert sample

        elif count == count_constant * 19:
            print("TRAJ 3 -----------  (LSPB)")

          
            pos_start = [0.0222, 1.4202, -1.0577, 0.0, 1.2955, 1.593] #B1 Approach
            pos_end   = [0.0222, 1.4565, -0.9375, 0.0, 1.1391, 1.593] #B1 Insert sample 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 11  —  B1 Insert sample → B1 Approach

        elif count == count_constant * 20:
            print("TRAJ 4 -----------  (LSPB )")

          
            pos_start = [0.0222, 1.4565, -0.9375, 0.0, 1.1391, 1.593] #B1 Insert sample
            pos_end   = [0.0222, 1.4202, -1.0577, 0.0, 1.2955, 1.593] #B1  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           
            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True



        # ------------------------------------------------------------------------
        
        elif count == count_constant * 21:
            print("Pos 9 -----------")
            position = [-0.5329, 1.4513, -1.0888, 0.0, 1.2956, -0.5329] #Lid holder approach
            piper_joint_movement(position)
            log_now = True

        # --------------------------------------------------------------------
        #  TRAJ 12  —   Lid holder approach → Lid holder Grip lid

        elif count == count_constant * 22:
            print("TRAJ 5 -----------  (LSPB )")

          
            pos_start = [-0.5329, 1.4513, -1.0888, 0.0, 1.2956, -0.5329] #Lid holder approach
            pos_end   = [-0.5329, 1.5664, -0.8417, 0.0, 0.9333, -0.5329] #Lid holder Grip lid


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 13  —   Lid holder approach → Lid holder Grip lid

        elif count == count_constant * 23:
            print("TRAJ 6 -----------  (LSPB )")

   
            pos_start = [-0.5329, 1.5664, -0.8417, 0.0, 0.9333, -0.5329] #Lid holder Grip lid
            pos_end   = [-0.5329, 1.4676, -1.0239, 0.0, 1.2144, -0.5329] #Lid holder approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            segment_durations = [3.0]  #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 24:
            print("Pos 12 -----------")
            position = [-0.5329, 1.2, -1.2, 0.0, 1.2144, -0.5329] # BO Incubator Lid on approach 
            piper_joint_movement(position)
            log_now = True

        elif count == count_constant * 25:
            print("Pos 12 -----------")
            position = [0.0111, 1.0171, -0.9242, -0.0046, 1.1725, 0.0088] # BO Incubator Lid on approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 14  —   BO Incubator Lid on approach  → BO Incubator Lid on

        elif count == count_constant * 26:
            print("TRAJ 7 -----------  (LSPB )")

           
            pos_start = [0.0111, 1.0171, -0.9242, -0.0046, 1.1725, 0.0088] #BO Incubator Lid on approach 
            pos_end   = [0.0075, 1.271, -0.8934, 0.0, 1.2805, 0.0075] #BO Incubator Lid on


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            # rad/s²

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        # --------------------------------------------------------------------
        #  TRAJ 15  —   BO Incubator Lid on  → BO Incubator Lid on approach

        elif count == count_constant * 27:
            print("TRAJ 8 -----------  (LSPB )")

          
            pos_start = [0.0075, 1.271, -0.8934, 0.0, 1.2805, 0.0075] #BO Incubator Lid on 
            pos_end   = [0.0073, 1.2739, -1.0104, 0.0000, 1.3945, 0.0073] #BO Incubator Lid on approach 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 28:
            print("Pos 15 -----------")
            position = [0.0, 0.9139, -0.6589, 0.0, -0.1678, 0.0] #Home position (Venter på at incubator er done med at do its thang)
            piper_joint_movement(position)
            log_now = True

        elif count == count_constant * 29:
            print("Pos 16 -----------")
            position = [0.0111, 1.0171, -0.9242, -0.0046, 1.1725, 0.0088] # BO Incubator Lid on approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 16  —   BO Incubator Lid on approach  → BO Incubator Lid on

        elif count == count_constant * 30:
            print("TRAJ 9 -----------  (LSPB )")

          
            pos_start = [0.0111, 1.0171, -0.9242, -0.0046, 1.1725, 0.0088] #BO Incubator Lid on approach 
            pos_end   = [0.0075, 1.271, -0.8934, 0.0, 1.2805, 0.0075] #BO Incubator Lid on


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        # --------------------------------------------------------------------
        #  TRAJ 17  —   BO Incubator Lid on  → BO Incubator Lid on approach

        elif count == count_constant * 31:
            print("TRAJ 10 -----------  (LSPB)")

           
            pos_start = [0.0075, 1.271, -0.8934, 0.0, 1.2805, 0.0075] #BO Incubator Lid on 
            pos_end   = [0.0111, 1.0171, -0.9242, -0.0046, 1.1725, 0.0088] #BO Incubator Lid on approach 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 32:
            print("Pos 19 -----------")
            position = [-0.5329, 1.4513, -1.0888, 0.0, 1.2956, -0.5329] #Lid holder approach
            piper_joint_movement(position)
            log_now = True

        # --------------------------------------------------------------------
        #  TRAJ 18 —   Lid holder approach → Lid holder Grip lid

        elif count == count_constant * 33:
            print("TRAJ 11 -----------  (LSPB )")

         
            pos_start = [-0.5329, 1.4513, -1.0888, 0.0, 1.2956, -0.5329] #Lid holder approach
            pos_end   = [-0.5226, 1.5575, -0.8406, 0.0, 0.9411, -0.5226] #Lid holder Grip lid


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 19  —   Lid holder approach → Lid holder Grip lid

        elif count == count_constant * 34:
            print("TRAJ 12 -----------  (LSPB )")

           
            pos_start = [-0.5226, 1.5575, -0.8406, 0.0, 0.9411, -0.5226] #Lid holder Grip lid
            pos_end   = [-0.5329, 1.4513, -1.0888, 0.0, 1.2956, -0.5329] #Lid holder approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

          
            segment_durations = [3.0]  #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        
        elif count == count_constant * 35:
            print("Pos 22 -----------")
            position = [0.0222, 1.4202, -1.0577, 0.0, 1.2955, 1.593] #B1 approach
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 20  —  B1 Approach → B1 Insert sample

        elif count == count_constant * 36:
            print("TRAJ 13 -----------  (LSPB )")

           
            pos_start = [0.0222, 1.4202, -1.0577, 0.0, 1.2955, 1.593] #B1 Approach
            pos_end   = [0.0222, 1.4623, -0.9244, 0.0, 1.1202, 1.593] #B1 Insert sample 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

         
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        # --------------------------------------------------------------------
        #  TRAJ 21 —  B1 Insert sample → B1 Approach

        elif count == count_constant * 37:
            print("TRAJ 14 -----------  (LSPB )")

            
            pos_start = [0.0222, 1.4623, -0.9244, 0.0, 1.1202, 1.593] #B1 Insert sample
            pos_end   = [0.0222, 1.4276, -1.0244, 0.0, 1.2548, 1.59] #B1  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

          
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        
        elif count == count_constant * 38:
            print("Pos 25 -----------")
            position = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray Approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 22  —  Sample tray approach → sample tray

        elif count == count_constant * 39:
            print("TRAJ 15 -----------  (LSPB)")

            
            pos_start = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach
            pos_end   = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 23  —  Sample tray→ sample tray approach 

        elif count == count_constant * 40:
            print("TRAJ 16 -----------  (LSPB)")

          
            pos_start = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray 
            pos_end   = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------


        elif count == count_constant * 41:
            print("Pos 22 -----------")
            position = [-0.1377, 1.2894, -0.8883, 0.0, 1.2569, -0.1526] #B5 approach
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 13  —  B5 Approach → B5 Insert sample

        elif count == count_constant * 42:
            print("TRAJ 24 -----------  (LSPB )")

            
            pos_start = [-0.1377, 1.2894, -0.8883, 0.0, 1.2569, -0.1526] #B5 Approach
            pos_end   = [-0.1526, 1.343, -0.8204, 0.0, 1.1354, -0.1526] #B5 Insert sample


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        # --------------------------------------------------------------------
        #  TRAJ 25  —  B5 Insert sample → B5 Approach

        elif count == count_constant * 43:
            print("TRAJ 14 -----------  (LSPB)")

            
            pos_start = [-0.1526, 1.343, -0.8204, 0.0, 1.1354, -0.1526] #B5 Insert sample
            pos_end   = [-0.1541, 1.2967, -0.93, 0.0, 1.3064, -0.1526] #B5  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        
        elif count == count_constant * 44:
            print("Pos 25 -----------")
            position = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray Approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 26  —  Sample tray approach → sample tray

        elif count == count_constant * 45:
            print("TRAJ 15 -----------  (LSPB)")

           
            pos_start = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach
            pos_end   = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray


            via_points = np.array([
                pos_start,
                pos_end,
            ])

            
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 27  —  Sample tray→ sample tray approach 

        elif count == count_constant * 46:
            print("TRAJ 16 -----------  (LSPB )")

        
            pos_start = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray 
            pos_end   = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

         
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 47:
            print("Pos 22 -----------")
            position = [0.1404, 1.2354, -0.845, 0.0, 1.2677, 2.356] #B10 approach
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 28  —  B10 Approach → B10 Insert sample

        elif count == count_constant * 48:
            print("TRAJ 13 -----------  (LSPB )")

        
            pos_start = [0.1404, 1.2354, -0.845, 0.0, 1.2677, 2.356] #B10 Approach
            pos_end   = [0.1315, 1.2696, -0.7583, 0.0, 1.1467, 2.356] #B10 Insert sample 


            via_points = np.array([
                pos_start,
                pos_end,
            ])

        
            segment_durations = [3.0]  #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        # --------------------------------------------------------------------
        #  TRAJ 29  —  B10 Insert sample → B10 Approach

        elif count == count_constant * 49:
            print("TRAJ 14 -----------  (LSPB )")

           
            pos_start = [0.1315, 1.2696, -0.7583, 0.0, 1.1467, 2.356] #B10 Insert sample
            pos_end   = [0.1315, 1.2359, -0.8677, 0.0, 1.2899, 2.356] #B10  Approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14            

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        
        elif count == count_constant * 50:
            print("Pos 25 -----------")
            position = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray Approach 
            piper_joint_movement(position)
            log_now = True
        # --------------------------------------------------------------------
        #  TRAJ 30  —  Sample tray approach → sample tray

        elif count == count_constant * 51:
            print("TRAJ 15 -----------  (LSPB )")

            
            pos_start = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach
            pos_end   = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray


            via_points = np.array([
                pos_start,
                pos_end,
            ])

          
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------
        # --------------------------------------------------------------------
        #  TRAJ 31  —  Sample tray→ sample tray approach 

        elif count == count_constant * 52:
            print("TRAJ 16 -----------  (LSPB )")

        
            pos_start = [-0.9601, 1.7088, -1.0975, 0.0, 1.0469, -0.9601] #Sample tray 
            pos_end   = [-0.9601, 1.6696, -1.1747, 0.0, 1.1632, -0.9601] #Sample tray approach


            via_points = np.array([
                pos_start,
                pos_end,
            ])

           
            segment_durations = [3.0]   #tid af bevægelsen fra A til B

            max_accel = 3.14           

            execute_lspb_trajectory(via_points, segment_durations, max_accel,
                                    dt=TRAJ_STREAM_DT)

            position = pos_end
            log_now = True
        # --------------------------------------------------------------------

        elif count == count_constant * 53:
            print("Pos 28 -----------")
            position = [0.0, 0.9139, -0.6589, 0.0, -0.1678, 0.0]   # Home Position
            piper_joint_movement(position)
            log_now = True

          
        elif count == count_constant * 54:
            print("Pos 29 -----------")
            position = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] #0 position (SLUT)
            piper_joint_movement(position)
            log_now = True


        if log_now:
            time.sleep(2)
            tcp_pose = piper.GetArmEndPoseMsgs()
            tcp_position = [
                tcp_pose.end_pose.X_axis,
                tcp_pose.end_pose.Y_axis,
                tcp_pose.end_pose.Z_axis,
                tcp_pose.end_pose.RX_axis,
                tcp_pose.end_pose.RY_axis,
                tcp_pose.end_pose.RZ_axis,
            ]
            log_to_csv(position, tcp_position)
            print(f"Logget nr {log_counter - 1}")
            time.sleep(1)

        time.sleep(0.005)
