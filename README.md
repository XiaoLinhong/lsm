### 前言

lsm: Lagrangian Stochastic Model

### 编译

lsm基于CMAKE进行编译，依赖CMAKE版本为3.x

```
# 第一步
mkdir build; cd build

# 第二步： 任选一个
FC=ifort cmake ../ # 普通编译
FC=ifort cmake -DCMAKE_BUILD_TYPE=Debug ../ # 调试

# 第三步
make
```

### 调试

```
gdb --args ./lsm.exe ../namelist.input
```

### 执行

```
./lsm.exe ../namelist.input

```
