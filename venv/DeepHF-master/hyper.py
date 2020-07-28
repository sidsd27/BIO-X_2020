import logging
import GPy
import GPyOpt
import matplotlib.pyplot as plt
plt.switch_backend('agg')

from training_util import *

#choose model to be trained
enzyme = 'esp'

# create logger with 'Model_application'
logger = logging.getLogger('Model')
logger.setLevel(logging.DEBUG)
# create file handler which logs even debug messages
fh = logging.FileHandler(enzyme + '_RNN_biofeat.log')
fh.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)

fc_activation_dict = {'1':'relu','2':'tanh', '3':'sigmoid', '4':'hard_sigmoid', '0':'elu'}
initializer_dict = {'1':'lecun_uniform','2':'normal', '3':'he_normal', '0':'he_uniform'}
optimizer_dict = {'1':SGD,'2':RMSprop, '3':Adagrad, '4':Adadelta,'5':Adam,'6':Adamax,'0':Nadam}
#parameter space
bounds = [
          #Discrete
          {'name':'em_drop','type':'discrete','domain':(0.1,0.2,0.3)},
          {'name':'rnn_drop','type':'discrete','domain':(0.1,0.2,0.4,0.4,0.5,0.6,0.7)},
          {'name':'rnn_rec_drop','type':'discrete','domain':(0.1,0.2,0.4,0.4,0.5,0.6,0.7)},
          {'name':'fc_drop','type':'discrete','domain':(0.1,0.2,0.4,0.4,0.5,0.6,0.7)},
          #Discrete
          {'name':'batch_size', 'type': 'discrete','domain': (20,30, 40, 50,60, 70, 80,90, 100)},
          {'name':'epochs', 'type': 'discrete','domain': (20, 25, 30, 35, 40, 45, 50, 55)},
          {'name':'em_dim', 'type': 'discrete','domain': (16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64, 68)},
          {'name':'rnn_units', 'type': 'discrete','domain': (20, 30, 40, 50, 60, 70, 80, 90, 100, 120, 140, 160, 180, 200, 220, 240)},
          {'name':'fc_num_hidden_layers', 'type': 'discrete','domain': (1,2,3,4,5,6)},
          {'name':'fc_num_units', 'type': 'discrete','domain': (50, 60, 70,80, 90, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360, 380, 400)},
          #Categorical
          {'name':'fc_activation', 'type': 'categorical','domain':tuple(map(int,tuple(fc_activation_dict.keys())))},
          {'name':'optimizer', 'type': 'categorical','domain':tuple(map(int,tuple(optimizer_dict.keys())))}
         ]

from dotmap import DotMap


def cat_decode(x):
    opt = DotMap()
    opt.em_drop = float(x[:, 0])
    opt.rnn_drop = float(x[:, 1])
    opt.rnn_rec_drop = float(x[:, 2])
    opt.fc_drop = float(x[:, 3])

    opt.batch_size = int(x[:, 4])
    opt.epochs = int(x[:, 5])
    opt.em_dim = int(x[:, 6])
    opt.rnn_units = int(x[:, 7])
    opt.fc_num_hidden_layers = int(x[:, 8])
    opt.fc_num_units = int(x[:, 9])

    opt.fc_activation = int(x[:, 10])
    opt.optimizer = int(x[:, 11])
    return opt
def f(x):
    print("inf")
    opt = cat_decode(x)
    param = {'model_type':enzyme,#mannualy set model_type as enzyme
            'em_drop':opt.em_drop,
           'rnn_drop':opt.rnn_drop,
           'rnn_rec_drop':opt.rnn_rec_drop,
           'fc_drop':opt.fc_drop,

           'batch_size':opt.batch_size,
           'epochs':opt.epochs,
           'em_dim':opt.em_dim,
           'rnn_units':opt.rnn_units,
           'fc_num_hidden_layers':opt.fc_num_hidden_layers,
           'fc_num_units':opt.fc_num_units,

           'fc_activation':opt.fc_activation,
           'optimizer':opt.optimizer}

    model = lstm_model(**param)
    y_test_pred = model.predict([X_test,X_test_biofeat])
    evaluation = mean_squared_error(y_test, y_test_pred)
    logger.info('----------')
    logger.info('params:{}'.format(param))
    logger.info('metrics——mse:{},spearman:{}'.format(evaluation,sp.stats.spearmanr(y_test, y_test_pred)[0]))
    return evaluation
opt_model = GPyOpt.methods.BayesianOptimization(f=f, domain=bounds, initial_design_numdata=30)
opt_model.run_optimization(max_iter = 300)
logger.info("optimized loss: {0}".format(opt_model.fx_opt))
for i,v in enumerate(bounds):
    name = v['name']
    logger.info('parameter {}:{}'.format(name,opt_model.x_opt[i]))
