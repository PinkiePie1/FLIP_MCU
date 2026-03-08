# 神秘嵌入式流体模拟器

受流体模拟项链项目的启发，用ch32v303复刻的跑在单片机上的物理模拟器。工匠精神手搓的代码到PIC仿真部分为止，后面计算压强和混合粒子速度的代码由codex生成。

## 目前功能

物理和算法部分已可在CH32V203上运转，摆脱了对FPU的依赖，纯模拟可跑到180FPS。firmware和对应的定点数代码详见[Here](https://github.com/PinkiePie1/CH32V203Dev/tree/master/Apps/led_fluid)。
该版本已接入LED矩阵和加速度传感器，成功实现了核心功能。

<img width="366" height="376" alt="image" src="https://github.com/user-attachments/assets/59a80178-057f-41ac-9ac1-86097e888443" />

当然，在有FPU的CH32v3上直接跑这个仓库里的单精度浮点数版本得快得多，由于CH32v303和V203几乎全部兼容，如需要进行迁移预计很轻松。

## 还需实现的部分
弄得更漂亮一点
