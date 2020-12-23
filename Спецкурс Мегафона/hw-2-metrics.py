import numpy as np

file = np.loadtxt('HW2_labels.txt', delimiter = ',')
y_predict, y_true = file[:,:2], file[:,-1]

def accuracy_score(y_true, y_predict, percent=None):
    
    if percent is None:
        p = 0.5   
        
    else:
        
        l = int(len(y_true)*percent/100.0)
        p = np.sort(y_predict[:,1])[-l]
    
    y_pred = (y_predict[:,1] >= p).astype(int)  
    y_true = y_true.astype(int) 
    result = np.mean(y_true == y_pred)
    
    return result

def precision_score(y_true, y_predict, percent=None):
    
    if percent is None:
        p = 0.5   
        
    else:
        
        l = int(len(y_true)*percent/100.0)
        p = np.sort(y_predict[:,1])[-l]
        
    
    y_pred = (y_predict[:,1] >= p).astype(int)  
    y_true = y_true.astype(int) 
    
    P = np.sum(y_pred)
    TP = sum([1 if y_pred[i] == y_true[i] & y_pred[i] == 1 else 0 for i in range(len(y_true))])

    result = 1.0*TP/P
    
    return result


def recall_score(y_true, y_predict, percent=None):
    
    if percent is None:
        p = 0.5   
        
    else:
        
        l = int(len(y_true)*percent/100.0)
        p = np.sort(y_predict[:,1])[-l]
        
    
    y_pred = (y_predict[:,1] >= p).astype(int)  
    y_true = y_true.astype(int) 
    
    P = np.sum(y_true)
    TP = sum([1 if y_pred[i] == y_true[i] & y_pred[i] == 1 else 0 for i in range(len(y_true))])

    result = 1.0*TP/P
    
    return result


def lift_score(y_true, y_predict, percent=None):
    
    if percent is None:
        p = 0.5   
        
    else:
        
        l = int(len(y_true)*percent/100.0)
        p = np.sort(y_predict[:,1])[-l]
        
    
    y_pred = (y_predict[:,1] >= p).astype(int)  
    y_true = y_true.astype(int)
    

    TP = sum([1 if y_pred[i] == y_true[i] & y_pred[i] == 1 else 0 for i in range(len(y_true))])

    result = 1.0*TP*len(y_true)/(np.sum(y_true)*np.sum(y_pred))
    
    return result   


def f1_score(y_true, y_predict, percent=None):
    
    if percent is None:
        p = 0.5   
        
    else:
        
        l = int(len(y_true)*percent/100.0)
        p = np.sort(y_predict[:,1])[-l]
        
    
    y_pred = (y_predict[:,1] >= p).astype(int)  
    y_true = y_true.astype(int) 
    

    TP = sum([1 if y_pred[i] == y_true[i] & y_pred[i] == 1 else 0 for i in range(len(y_true))])

    result = 2.0*TP/(np.sum(y_true)+np.sum(y_pred))
    
    return result   


f1_score(y_true, y_predict)              