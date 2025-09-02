import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import re
import BZdrawer
plt.rcParams['font.family'] = 'Times New Roman'
# important parameter declare
vb = 27
cb = 28


#7.2pbe
E_VB = 13.0210

b1 = np.array([ 1.000000,  0.577350,  -0.000000])
b2 = np.array([ 0.000000,  1.154701,  0.000000])
b3 = np.array([ 0.000000,  -0.000000,  0.620690])


nk_x=20
nk_y=20
nk_z=20

kvectors = np.vstack((b1, b2, b3))
transformation_matrix = np.vstack((b1, b2, b3)).T

# high symmetry point
Ψ = 0.73365
ζ = 0.39699
η = 0.58951
φ = 0.74186
high_sym_point = {
    "Γ": np.array([0.0, 0.0, 0.0]),
    "M": np.array([0.5, 0.0, 0.0]),
    "K": np.array([0.333, 0.333, 0.0]),
    "A": np.array([0.0, 0.0, 0.5]),
    "L": np.array([0.5, 0.0, 0.5]),
    "H": np.array([0.333, 0.333, 0.5]),
}



xdiv=1/nk_x
ydiv=1/nk_y
zdiv=1/nk_z
# new : 找valley
def found_valley(crys_coord,energy):

    valleys = {}
    unit_crys_coord_dict = {}
    crys_coord_dict = {}

    for i in range(len(crys_coord)):
        unit_crys_coord_dict[crys_coord[i]] = energy[i]
        x = crys_coord[i][0]
        y = crys_coord[i][1]
        z = crys_coord[i][2]
        unit_crys_coord_dict[(-x,-y,-z)] = energy[i]
        
    print(f'unit_crys_coord_dict:{len(unit_crys_coord_dict)}')

    cut_x_max = 1
    cut_x_min = -1
    cut_y_max = 1
    cut_y_min = -1
    cut_z_max = 1
    cut_z_min = -1
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for coor, energy in unit_crys_coord_dict.items(): 
                    x = coor[0] + (-1 + i)
                    y = coor[1] + (-1 + j)
                    z = coor[2] + (-1 + k)

                    if (x <= cut_x_max and x >= cut_x_min) \
                        and (y <= cut_y_max and y >= cut_y_min) \
                        and (z <= cut_z_max and z >= cut_z_min):
                        crys_coord_dict[(x,y,z)] = energy

                        # 轉座標，畫圖用
                        # coordinate = np.array([x,y,z])
                        # coordinate = np.dot(transformation_matrix, coordinate)
                        # coordinate_x.append(coordinate[0])
                        # coordinate_y.append(coordinate[1])
                        # coordinate_z.append(coordinate[2])
                        # enengy_crys.append(energy)


    print(f'crys_coord_dict:{len(crys_coord_dict)}')

    # 取範圍
    x_range = int((cut_x_max - cut_x_min)/xdiv + 1)
    y_range = int((cut_y_max - cut_y_min)/ydiv + 1)
    z_range = int((cut_z_max - cut_z_min)/zdiv + 1)

    # 創空三維陣列
    array_shape = (x_range, y_range, z_range)
    xyz_array = np.zeros(array_shape)    
    print("xyz_array shape:", xyz_array.shape)

    # 把字典中的值填到三維陣列中的對應位置
    for xyz_coord, value in crys_coord_dict.items():

        # 不round會出錯
        x_idx = int(round((xyz_coord[0] - cut_x_min) / xdiv, 1))  # 計算 x 座標的索引
        y_idx = int(round((xyz_coord[1] - cut_y_min) / ydiv, 1))  # 計算 y 座標的索引
        z_idx = int(round((xyz_coord[2] - cut_z_min) / zdiv, 1))  # 計算 z 座標的索引
        xyz_array[x_idx, y_idx, z_idx] = value


    # 從陣列抓出valley
    for z in range(z_range-2):
        for y in range(y_range-2):
            for x in range(x_range-2):
                # 中心位置(center)
                x_c = x + 1
                y_c = y + 1
                z_c = z + 1

                # 從3x3x3確認中心是否最小
                min_check_flag = True
                for i in [-1,0,1]:
                    for j in [-1,0,1]:
                        for k in [-1,0,1]:
                            # 周圍(surround)
                            x_s = x_c + i
                            y_s = y_c + j
                            z_s = z_c + k
                            if (i == 0 and j == 0 and k == 0):
                                continue
                            else:
                                ## 若比周遭任何一個大就False，找VB valley用>= ，找CB local max.用<=
                                if (xyz_array[x_c,y_c,z_c] >= xyz_array[x_s,y_s,z_s]):
                                    min_check_flag = False
                                # if ([x_c,y_c,z_c] == [20,20,20]):
                                #     print(x_s,y_s,z_s,":",xyz_array[x_s,y_s,z_s] )


                # 新增valleys
                if (min_check_flag == True):
                    if (xyz_array[x_c,y_c,z_c] < 5):
                        crys_x_c = x_c*xdiv + cut_x_min
                        crys_y_c = y_c*ydiv + cut_x_min
                        crys_z_c = z_c*zdiv + cut_x_min
                        
                        valleys[(crys_x_c,crys_y_c,crys_z_c)] = round(xyz_array[x_c,y_c,z_c],4)
    
    print(f'valleys={valleys}')
    print(f'center={xyz_array[nk_x,nk_y,nk_z]}')

    Cart_valleys = {}
        
    # x : 1.3, y : 2.2, z : 1.3
    cut_x_max = b2[1]/np.sqrt(3)*2 *0.5 +0.1
    cut_x_min = -cut_x_max
    cut_y_max = b2[1]*0.5+0.1 
    cut_y_min = -b2[1]*0.5-0.1
    cut_z_max = b3[2]*0.5
    cut_z_min = -b3[2]*0.5

    for key in valleys.keys():
        crys_coordinate = np.array([key[0],key[1],key[2]])
        Cart_coordinate = np.dot(transformation_matrix, crys_coordinate)
        x,y,z = Cart_coordinate
        # 把BZ以外的刪掉
        if (x <= cut_x_max and x >= cut_x_min) \
            and (y <= cut_y_max and y >= cut_y_min) \
            and (z <= cut_z_max and z >= cut_z_min):
            Cart_valleys[(x,y,z)] = valleys[key]
        
    print(f'Cart_valleys={Cart_valleys}')
    return Cart_valleys


##目前沒用到
def isosurface(Cart_x, Cart_y, Cart_z, energys):
    isosurface_point = {}
    surface_E_max = 6.0
    surface_E_min = 5.8
    for i,energy in enumerate(energys):
        if (energy >= surface_E_min and energy <= surface_E_max):
            isosurface_point[(Cart_x[i],Cart_y[i],Cart_z[i])] = energy
    return isosurface_point




count = 0

# 讀取資料
f = open("./nscf_ZnO.out", mode = 'r')
read_energy_flag = False
coordinates_data = {}
cryst_coord_flag = False 
cryst_coord = []
cc = 0
for line in f:
    if line.strip() == "cryst. coord.":
        cryst_coord_flag = True 

    match = re.match(r'k\(\s*\d+\s*\) = \(\s*([-+]?\d+\.\d+)\s*([-+]?\d+\.\d+)\s*([-+]?\d+\.\d+)\s*\)', line.strip())
    if cryst_coord_flag == True :
        if match:
            x, y, z = match.groups()
            cryst_coord.append((float(x), float(y), float(z)))
            cc = cc + 1

    match = re.match(r'k =\s*([-+]?[0-9]+\.[0-9]{4})\s*([-+]?[0-9]+\.[0-9]{4})\s*([-+]?[0-9]+\.[0-9]{4})', line.strip())
    if match:
        # 將match到的三位小數取出
        x, y, z = match.groups()
        count = count + 1
        # 做成tuple再存入字典當作key,因key不可用list這種會變的
        current_coordinate = (float(x), float(y), float(z))
        coordinates_data[current_coordinate] = []
 
        read_energy_flag = True #開啟偵測能量
        continue

    if read_energy_flag == True:
        if line.strip() == "": #跳過空的
            continue
        elif line.strip() == "occupation numbers": #看到這個就關掉
            read_energy_flag = False
        else:
            numbers_str = re.findall(r'[-+]?\d+\.\d+', line.strip())
            if numbers_str:
                numbers = list(map(float, numbers_str))
                coordinates_data[current_coordinate].extend(numbers)


# 用正則抓所有浮點數，避免數字黏在一起
numbers_str = re.findall(r'[-+]?\d+\.\d+', line.strip())
if numbers_str:  # 確保有抓到
    numbers = list(map(float, numbers_str))
    coordinates_data[current_coordinate].extend(numbers)

# print(count)
# print(cc)


coordinate_x = []
coordinate_y = []
coordinate_z = []
energy_VB = []
energy_CB = []
coordinates = []

for coordinate, energy in coordinates_data.items():

    # coordinate = np.dot(np.linalg.inv(transformation_matrix), coordinate)
    coordinate_x.append(coordinate[0])
    coordinate_y.append(coordinate[1])
    coordinate_z.append(coordinate[2])
    coordinates.append(coordinate)
    VB = energy[vb-1] - E_VB
    CB = energy[cb-1] - E_VB
    energy_VB.append(VB)
    energy_CB.append(CB)

# print(min_Eg)

## 找valley
# input : relative coordinate, arbitary enengy, output : valleys{keys = Cart_coordinate : value = energy}
valleys = found_valley(cryst_coord,energy_CB)
valleys_x = [i[0] for i in valleys.keys()]
valleys_y = [i[1] for i in valleys.keys()]
valleys_z = [i[2] for i in valleys.keys()]

# 擴充k point
coordinate_x = np.array(coordinate_x)
coordinate_y = np.array(coordinate_y)
coordinate_z = np.array(coordinate_z)
coordinates = np.array(coordinates)

unit_coordinate_x = np.concatenate([coordinate_x,-coordinate_x])
unit_coordinate_y = np.concatenate([coordinate_y,-coordinate_y])
unit_coordinate_z = np.concatenate([coordinate_z,-coordinate_z])
unit_coordinates = np.concatenate([coordinates,-coordinates])
energy_VB = energy_VB*2
energy_CB = energy_CB*2


total_coordinate_x = []
total_coordinate_y = []
total_coordinate_z = []
total_energy_VB = []
total_energy_CB = []


# x : 1.3, y : 2.2, z : 1.15
# cut_x_max = 0.2
# cut_x_min = -1.2
# cut_y_max = 2.7
# cut_y_min = 1.4
# cut_z_max = 1.0
# cut_z_min = -1.0
#下面是turn =1 是學長設定

turn = 0  # turn 1 有完整剖面  turn 0 自己設定只看z 點 因G跟A (x,y,z)的 x,y重疊


if(turn ==1 ) :
  cut_x_max = b2[1]/np.sqrt(3)*2 *0.5 +0.1
  cut_x_min = -cut_x_max
  cut_y_max = b2[1]*0.5+0.1 
  cut_y_min = -b2[1]*0.5-0.1
  cut_z_max = 0.05
  cut_z_min = -0.05
if(turn == 0 ) :
  epsilon = 0.05   # 容許 X,Y 的小範圍
  cut_x_max = epsilon
  cut_x_min = -epsilon
  cut_y_max = epsilon
  cut_y_min = -epsilon
  cut_z_max = b3[2] * 0.5 + 1e-3   # 比 A 再大一點，避免切掉
  cut_z_min = -0.1          # 可以多留一點 Γ 附近





ak = 2   # 1 = 完整剖面, 0 = 只看 Γ–A

# ---------- 先設定初始 cut ----------
if ak == 1:
    # 想要看的大範圍
    cut_x_max = b2[1]/np.sqrt(3)*2 *0.5 + 0.1
    cut_x_min = -cut_x_max
    cut_y_max = b2[1]*0.5 + 0.5
    cut_y_min = -b2[1]*0.5 - 0.5
    cut_z_max = b3[2]*0.5 + 0.5
    cut_z_min = -b3[2]*0.5 - 0.5
    bz_x = b2[1]/np.sqrt(3)*0.5 # ---------- 再切回布里淵區 ----------
    bz_y = b2[1]*0.5  # ---------- 再切回布里淵區 ----------
    bz_z = b3[2]*0.5# ---------- 再切回布里淵區 ----------
    
    cut_x_max = min(cut_x_max,  bz_x)
    cut_x_min = max(cut_x_min, -bz_x)
    cut_y_max = min(cut_y_max,  bz_y)
    cut_y_min = max(cut_y_min, -bz_y)
    cut_z_max = min(cut_z_max,  bz_z)
    cut_z_min = max(cut_z_min, -bz_z)

if ak == 0:
    epsilon = 0.05
    cut_x_max = epsilon
    cut_x_min = -epsilon
    cut_y_max = epsilon
    cut_y_min = -epsilon
    cut_z_max = b3[2]*0.5 + 0.1   # 先放大一點
    cut_z_min = -0.1              # 包含 Γ
    
    cut_x_max = min(cut_x_max,  bz_x)
    cut_x_min = max(cut_x_min, -bz_x)
    cut_y_max = min(cut_y_max,  bz_y)
    cut_y_min = max(cut_y_min, -bz_y)
    cut_z_max = min(cut_z_max,  bz_z)
    cut_z_min = max(cut_z_min, -bz_z)

# ---------- 再切回布里淵區 ----------











for i in range(3):
    for j in range(3):
        for k in range(3):
            for index, c in enumerate(unit_coordinates):
                shift = np.array(c) + b1 * (-1 + i) + b2 * (-1 + j) + b3 * (-1 + k) 
                if (shift[0] <= cut_x_max and shift[0] >= cut_x_min) \
                      and (shift[1] <= cut_y_max and shift[1] >= cut_y_min) \
                      and (shift[2] <= cut_z_max and shift[2] >= cut_z_min):
                    total_coordinate_x.append(shift[0])
                    total_coordinate_y.append(shift[1])
                    total_coordinate_z.append(shift[2])
                    total_energy_VB.append(energy_VB[index])
                    total_energy_CB.append(energy_CB[index])

isosurface_point = isosurface(total_coordinate_x,total_coordinate_y,total_coordinate_z,total_energy_CB)
isosurface_point_x = [i[0] for i in isosurface_point.keys()]
isosurface_point_y = [i[1] for i in isosurface_point.keys()]
isosurface_point_z = [i[2] for i in isosurface_point.keys()]
isosurface_energy = [i for i in isosurface_point.values()]


## 畫圖相關
fig = plt.figure()
fig2 = plt.figure()

ax = fig.add_subplot(111, projection='3d')
bx = fig2.add_subplot(111, projection='3d')
BZ1 = BZdrawer.BZ(kvectors)
BZ1.bulkBZ()
BZ1.draw_bulkBZ(ax,bx)
# sc = ax.scatter(unit_coordinate_x, unit_coordinate_y, unit_coordinate_z, c=energy_CB, cmap='terrain', marker='o', alpha=0.2)
# sc = ax.scatter(coordinate_x_crys,coordinate_y_crys,coordinate_z_crys,c = enengy_crys, cmap='terrain', marker='o', alpha=0.2)
sc = ax.scatter(total_coordinate_x, total_coordinate_y, total_coordinate_z, c=total_energy_CB, cmap='gnuplot', marker='o', alpha=1) #alpha為透明度
# sc = ax.scatter(isosurface_point_x, isosurface_point_y, isosurface_point_z, c=isosurface_energy, cmap='gnuplot')

# 畫valley位置(原本)
bx.scatter(valleys_x, valleys_y, valleys_z,c='purple', s=100, marker='x', label='Valleys')
#cl=1 標示出valley座標
cl=1
if(cl==1):
# 在 bx 標出座標與能量
 for (x, y, z), E in valleys.items():
    bx.text(x, y, z+0.03,  # +0.03 是避免文字壓在點上
            f"({x:.2f},{y:.2f},{z:.2f})\nE={E:.2f} eV",
            fontsize=8, color='purple', ha='center')


# 畫 valley 位置並標上座標
# 畫 valley 位置 (右邊圖 bx)
# 畫 valley 位置 (右邊圖 ax)
ax.scatter(valleys_x, valleys_y, valleys_z,
           c='purple', s=100, marker='x', label='Valleys')

# 在右邊圖 ax 標出座標與能量 choode 是看3d座標點
choose =0 
if(choose==1) : 
 for (x, y, z), E in valleys.items():
    ax.text(x, y, z+0.05,
            f"({x:.2f},{y:.2f},{z:.2f})\nE={E:.2f} eV",
            fontsize=10, color='purple', ha='center')





# ax.scatter(0, 0, 0, c='k', marker='o', label='Origin')
ax.quiver(0, 0, 0, b1[0], b1[1], b1[2], color='r', label='b1', arrow_length_ratio=0.1)
ax.quiver(0, 0, 0, b2[0], b2[1], b2[2], color='g', label='b2', arrow_length_ratio=0.1)
ax.quiver(0, 0, 0, b3[0], b3[1], b3[2], color='b', label='b3', arrow_length_ratio=0.1)
# bx.scatter(0, 0, 0, c='k', marker='o', label='Origin')
bx.quiver(0, 0, 0, b1[0], b1[1], b1[2], color='r', label='b1', arrow_length_ratio=0.1)
bx.quiver(0, 0, 0, b2[0], b2[1], b2[2], color='g', label='b2', arrow_length_ratio=0.1)
bx.quiver(0, 0, 0, b3[0], b3[1], b3[2], color='b', label='b3', arrow_length_ratio=0.1)
ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')
bx.set_xlabel('X Label')
bx.set_ylabel('Y Label')
bx.set_zlabel('Z Label')

cart_sym_point_x = []
cart_sym_point_y = []
cart_sym_point_z = []
# high symmetry point from fraction convert to Cart.
for name, frac_sym_point in high_sym_point.items():
    cart_sym_point = np.dot(frac_sym_point, kvectors)
    cart_sym_point_x.append(cart_sym_point[0])
    cart_sym_point_y.append(cart_sym_point[1])
    cart_sym_point_z.append(cart_sym_point[2])
    ax.text(cart_sym_point[0], cart_sym_point[1], cart_sym_point[2], name, fontsize=20)
    ax.scatter(cart_sym_point[0], cart_sym_point[1], cart_sym_point[2], c='k', marker='o')
    bx.text(cart_sym_point[0], cart_sym_point[1], cart_sym_point[2], name, fontsize=20, zorder=10)
    bx.scatter(cart_sym_point[0], cart_sym_point[1], cart_sym_point[2], c='k', marker='o')

# 隱藏坐標軸 (axis lines, ticks, and labels)
for axes in [bx]:
    axes.set_axis_off()  # 通用方式，完全隱藏整個 3D 坐標系
#改觀看角度

bx.view_init(elev=60, azim=-20) 
# print(BZ1.hs_points)
cbar = plt.colorbar(sc, label='energy_CB')
ax.legend(loc='upper left', fontsize='16')
bx.legend(loc='upper left', fontsize='16')

axe_min = -0.7
axe_max = 0.7
ax.set_xlim(axe_min, axe_max)
ax.set_ylim(axe_min, axe_max)
ax.set_zlim(axe_min, axe_max)
bx.set_xlim(axe_min, axe_max)
bx.set_ylim(axe_min, axe_max)
bx.set_zlim(axe_min, axe_max)
sc.set_clim(min(energy_CB), max(energy_CB))
plt.tight_layout()


plt.show()
# plt.savefig('output.png')


    