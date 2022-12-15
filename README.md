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
Sample Result:
![image](https://user-images.githubusercontent.com/85045910/207783785-95665491-7642-4f7c-8aa1-cd8deafb4b69.png)
![image](https://user-images.githubusercontent.com/85045910/207783920-66b06b1f-394f-40e9-83c3-7f541c2b010b.png)
![image](https://user-images.githubusercontent.com/85045910/207784084-dcd4292a-b5ef-430b-bfa8-abe4f2b22ee8.png)

