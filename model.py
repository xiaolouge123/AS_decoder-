# coding: utf-8

import numpy as np
from keras.preprocessing.text import Tokenizer, one_hot
from keras.preprocessing.sequence import pad_sequences
from keras.utils import to_categorical
from keras.layers import Dense, Input, GlobalMaxPooling1D
from keras.layers import Conv1D, MaxPooling1D, Embedding, LSTM, Flatten
from keras.models import Model


GLOVE_DIR =  './data/nucle.txt'
TEXT_DATA_DIR = './data/data_to_train.txt'
MAX_SEQUENCE_LENGTH = 120
MAX_NB_WORDS = 20
EMBEDDING_DIM = 5
VALIDATION_SPLIT = 0.8
CATEGORY = {'AS junction' : 1,
            'non junction' : 2}

# first, build index mapping words in the embeddings set to their embedding vector
print('Indexing word vectors.')

embeddings_index = {}
f = open(GLOVE_DIR, 'r')
for line in f:
    values = line.strip().split(' ')
    word = values[0]
    coefs = np.asarray(values[1:], dtype='float32')
    embeddings_index[word] = coefs
f.close()

print('Found %s word vectors.' % len(embeddings_index))

# second, prepare text samples and their labels
print('Processing text dataset')

texts = []  # list of text samples
labels_index = CATEGORY  # dictionary mapping label name to numeric id
labels = []  # list of label ids
f = open(TEXT_DATA_DIR, 'r', encoding='utf8')
lines = f.readlines()
f.close()
for line in lines:
    tmp = line.strip().split('\t')
    AS_data = ''
    for seq in tmp[3:7]:
        AS_data += seq
    texts.append(AS_data)
    labels.append(int(tmp[-1]))

print('Found %s texts.' % len(texts))

# finally, vectorize the text samples into a 2D integer tensor
tokenizer = Tokenizer(num_words=MAX_NB_WORDS,char_level=True)
tokenizer.fit_on_texts(texts)
sequences = tokenizer.texts_to_sequences(texts)
print(sequences[0])
word_index = tokenizer.word_index
print('Found %s unique tokens.' % len(word_index))
for i in word_index:
    print(i,word_index[i])

data = pad_sequences(sequences, maxlen=MAX_SEQUENCE_LENGTH,value=1)
print(data[0])

labels = to_categorical(np.asarray(labels))
print('Shape of data tensor:', data.shape)
print('Shape of label tensor:', labels.shape)

# split the data into a training set and a validation set
indices = np.arange(data.shape[0])
np.random.shuffle(indices)
data = data[indices]
labels = labels[indices]
num_validation_samples = int(VALIDATION_SPLIT * data.shape[0])

x_train = data[:num_validation_samples]
y_train = labels[:num_validation_samples]
x_val = data[num_validation_samples:]
y_val = labels[num_validation_samples:]


# prepare embedding matrix
print('Preparing embedding matrix.')
num_words = min(MAX_NB_WORDS, len(word_index)) + 1
embedding_matrix = np.zeros((num_words, EMBEDDING_DIM))
error = 0
random_error = 0
for word, i in word_index.items():
    # print(word,i)
    if word not in embeddings_index:
        error += 1
    if i >= MAX_NB_WORDS:
        continue
    embedding_vector = embeddings_index.get(word)
    if embedding_vector is not None:
        # words not found in embedding index will be all-zeros.
        embedding_matrix[i] = embedding_vector
    else :
        embedding_matrix[i] = np.random.uniform(-0.01, 0.01, 5)
        random_error += 1
print('Out of vocabulary: ', error)
print('Random embedding: ', random_error)

# load pre-trained word embeddings into an Embedding layer
# note that we set trainable = False so as to keep the embeddings fixed
embedding_layer = Embedding(num_words,
                            EMBEDDING_DIM,
                            weights=[embedding_matrix],
                            input_length=MAX_SEQUENCE_LENGTH,
                            trainable=False)


print('Training model.')

# train a 1D convnet with global maxpooling
sequence_input = Input(shape=(MAX_SEQUENCE_LENGTH,), dtype='float64')
embedded_sequences = embedding_layer(sequence_input)
# x = LSTM(128,dropout=0.2,recurrent_dropout=0.2)(embedded_sequences)
x = Conv1D(128, 5, activation='relu')(embedded_sequences)
x = MaxPooling1D(5)(x)
# x = Conv1D(128, 5, activation='relu')(x)
# x = MaxPooling1D(5)(x)
x = Conv1D(128, 5, activation='relu')(x)
x = GlobalMaxPooling1D()(x)
x = Dense(128, activation='relu')(x)
preds = Dense(len(labels_index), activation='sigmoid')(x)

model = Model(sequence_input, preds)
model.compile(loss='binary_crossentropy',
              optimizer='rmsprop',
              metrics=['acc'])

model.fit(x_train, y_train,
          batch_size=512,
          epochs=10,
          validation_data=(x_val, y_val))

word_index_reverse = {}
for i in word_index:
    word_index_reverse[word_index[i]] = i


r = model.predict(x_val)
pred_result = []
for item in list(r):
    max_i = 0
    max_index = 0
    for index, i in enumerate(item):
        if i > max_i:
            max_i = i
            max_index = index
    pred_result.append(max_index + 1)
pred_true = list(int(x) for x in y_val.dot(np.array([1, 0]).T))
f = open('./data/results.txt', 'w', encoding='utf8')
assert len(pred_result) == len(pred_true)
for i in range(len(pred_true)):
    s = ''
    for w in list(x_val)[i]:
        if w != 0:
            if word_index_reverse.get(w) != None:
                s += word_index_reverse[w]
    s += '\t' + str(pred_true[i]) + '\t' + str(pred_result[i]) + '\n'
    f.write(s)
f.close()

 