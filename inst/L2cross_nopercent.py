def L2cross_nopercent(path, coxnnet_method, percent, nfold):
	import numpy
	import sklearn
	import theano
	from sklearn.model_selection import train_test_split

	import sys
	import os
	
	sys.path.append(path)
	from cox_nnet3 import cox_nnet3

	# load data
	os.chdir(path)
	x = numpy.loadtxt(fname="x.csv",delimiter=",",skiprows=0)
	ytime = numpy.loadtxt(fname="ytime.csv",delimiter=",",skiprows=0)
	ystatus = numpy.loadtxt(fname="ystatus.csv",delimiter=",",skiprows=0)

	# split into test/train sets
	# x_train, x_test, ytime_train, ytime_test, ystatus_train, ystatus_test = \
	# 	train_test_split(x, ytime, ystatus, test_size = 1-percent, random_state = 100)

	#Define parameters
	model_params = dict(node_map = None, input_split = None)
	search_params = dict(method = coxnnet_method, learning_rate=0.01, momentum=0.9,
	    max_iter=2000, stop_threshold=0.995, patience=1000, patience_incr=2, rand_seed = 123,
	    eval_step=23, lr_decay = 0.9, lr_growth = 1.0)
	cv_params = dict(cv_seed=1, n_folds=int(nfold), cv_metric = "loglikelihood", L2_range = numpy.arange(-4.5,1,0.5))

	#cross validate training set to determine lambda parameters
	cv_likelihoods, L2_reg_params, mean_cvpl = cox_nnet3.L2CVProfile(x,ytime,ystatus,
	    model_params,search_params,cv_params, verbose=False)

	numpy.savetxt("ttcv_likelihoods.csv", cv_likelihoods, delimiter=",")


	#build model based on optimal lambda parameter
	L2_reg = L2_reg_params[numpy.argmax(mean_cvpl)]
	model_params = dict(node_map = None, input_split = None, L2_reg=numpy.exp(L2_reg))
	model, cost_iter = cox_nnet3.trainCoxMlp(x, ytime, ystatus, model_params, search_params, verbose=True)

	theta = model.predictNewData(x)

	numpy.savetxt("theta.csv", theta, delimiter=",")
	numpy.savetxt("ytime.csv", ytime, delimiter=",")
	numpy.savetxt("ystatus.csv", ystatus, delimiter=",")

	
