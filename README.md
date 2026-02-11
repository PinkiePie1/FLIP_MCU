# 神秘嵌入式流体模拟器
受流体模拟项链项目的启发，用ch32v303复刻的跑在单片机上的物理模拟器

## 目前功能

物理和算法部分已可在CH32V203上运转，摆脱了对FPU的依赖，纯模拟可跑到130FPS。代码详见[此处](https://github.com/PinkiePie1/ch32DevEnv/tree/master/CH32V203apps/FLIP)：

![test](https://github.com/user-attachments/assets/98db300b-1680-49b1-b5f7-09cf8e177358)

当然，在有FPU的CH32v3上跑得要快得多。

## 还需实现的部分

把输入改为重力感应
显示更改为LED点阵，基于DMA的LED点阵驱动。
GPIO够多，可用共阳点阵屏
