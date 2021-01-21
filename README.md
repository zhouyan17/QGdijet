# qgdijet
Roo.C用来从workspace的结果中提取histogram，Altasstyle的部分是为了画图检验的，没什么用可以删了。
直接root Roo.C就可以

bumphunter已经编译好了，直接跑这个就行：
./runBumpHunter --inFile Roo.root --outPath results/ --outFileName a.root --dataHist basicData --bkgHist basicBkg --nPseudoExpBH 1000 --minBH 1000 --maxBH 10000
--inFile Roo.root --outPath results/ --outFileName a.root用来控制input和output的名字

bumphunter的结果会存在root里，再跑root draw.C画图，需要把bumphunter输出在屏幕上的结果作为draw.C的参数，相应的total和GG情况的参数我写在最后的注释里了。
