# Cache Placement Reproduce
## 簡介
重現 [Cache Placement in Two-Tier HetNets With Limited Storage Capacity: Cache or Buffer?](https://ieeexplore.ieee.org/document/8382257) 這篇文獻的實驗
## 進度
1. 完成求Smf會需要的參數和方程式
2. (WIP)迭代Smf找最適合的值
### 目前遇到的問題
1. Smf的數值不合理（小於0且趨近於0）
    1.1 已檢查係數正確
    1.2 (WIP)檢查S_root是否合理 (目前小於0)
