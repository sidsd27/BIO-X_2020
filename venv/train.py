from prediction_util import *
from training_util import *
param = {'model_type':'esp',
        'em_drop':0.2,
        'rnn_drop':0.5,
        'rnn_rec_drop':0.4,
        'fc_drop':0.4,
        'batch_size':80,
        'epochs':45,
        'em_dim':44,
        'rnn_units':80,
        'fc_num_hidden_layers':3,
        'fc_num_units':300,
        'fc_activation':'3',
        'optimizer':'6'}

trained_esp_model = lstm_model(**param)
get_metrics(trained_esp_model)
param = {'model_type':'hf',
    'em_drop': 0.2, 
     'rnn_drop': 0.5, 
     'rnn_rec_drop': 0.4, 
     'fc_drop': 0.5, 
     'batch_size': 80, 
     'epochs': 45, 
     'em_dim': 48, 
     'rnn_units': 80, 
     'fc_num_hidden_layers': 2, 
     'fc_num_units': 300, 
     'fc_activation': 3, 
     'optimizer': 6}

trained_hf_model = lstm_model(**param)
get_metrics(trained_hf_model,'hf')
