# HR_diagram_generator
INPUTS:
1. B band image (.fits)
2. V band image (.fits)
3. Star list data fetch from DS9 (.csv)
4. Ground truth data from  Stellar Evolution Database
Astronomy_3.m is the Final code
已做更新
1.  利用三圓圈來確認星點數量，是否一次只偵測到一顆星，且利用下圖方式計算counts of a star，reference from Prof. Shih-Ping Lai's lecture
![image](https://user-images.githubusercontent.com/85045910/148628209-0cf3ca0c-f557-4fff-b37c-20a2f1023583.png)
2. 偵測到最大亮度會告知且不取用此星點

