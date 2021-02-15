challengeutils query "select * from evaluation_9614582" > coarse_post.csv
challengeutils query "select * from evaluation_9614583" --limit 10 > fine_post1.csv
challengeutils query "select * from evaluation_9614583" --limit 10 --offset 10 > fine_post2.csv
challengeutils query "select * from evaluation_9614583" --limit 10 --offset 20 > fine_post3.csv
challengeutils query "select * from evaluation_9614583" --limit 10 --offset 30 > fine_post4.csv
