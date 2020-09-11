challengeutils query "select * from evaluation_9614317" > coarse_final.csv
challengeutils query "select * from evaluation_9614318" --limit 10 > fine_final1.csv
challengeutils query "select * from evaluation_9614318" --limit 10 --offset 10 > fine_final2.csv
challengeutils query "select * from evaluation_9614318" --limit 10 --offset 20 > fine_final3.csv
challengeutils query "select * from evaluation_9614318" --limit 10 --offset 30 > fine_final4.csv
challengeutils query "select * from evaluation_9614318" --limit 10 --offset 40 > fine_final5.csv
