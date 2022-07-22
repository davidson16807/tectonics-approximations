import csv
import numpy as np
import tensorflow as tf
import plotly.express as px

with open('varsanyi1986.csv', newline='') as csvfile:
    csvreader = csv.reader(csvfile, delimiter=',', quotechar='"')
    next(csvreader, None) # skip header
    observations = [[float(Zc),float(p),float(T), float(Z)]
                    for (i,p,Zc,T,Z) in csvreader
                    ]
                  # if p >= 5 or (p>=1.2 and T<1.05) ]

I = np.array([[Zc, p, 1/T, p/T] for Zc,p,T,Z in observations])
V = np.array([[p/T] for Zc,p,T,Z in observations])
Tinv = np.array([[1/T] for Zc,p,T,Z in observations])
Zc = np.array([[Zc] for Zc,p,T,Z in observations])
Z = np.array([[Z] for Zc,p,T,Z in observations])

dropout_change_per_layer = 0.02
model = tf.keras.Sequential([
    # tf.keras.layers.BiasLayer(),
    # tf.keras.layers.Dropout(dropout_change_per_layer*1),
    # tf.keras.layers.BiasLayer(),
    tf.keras.layers.Dense(4, activation='relu'),
    tf.keras.layers.Dense(4, activation='sigmoid'),
    # tf.keras.layers.Dropout(dropout_change_per_layer*1),
    # tf.keras.layers.Dropout(dropout_change_per_layer*0),
    # tf.keras.layers.BiasLayer(),
    # tf.keras.layers.Dense(2, activation='relu'),
    # tf.keras.layers.Dropout(dropout_change_per_layer*0),
    # tf.keras.layers.BiasLayer(),
    tf.keras.layers.Dense(1, activation='linear'),
    # tf.keras.layers.Dropout(dropout_change_per_layer*0),
    # tf.keras.layers.BiasLayer(),
])
model.compile(optimizer=tf.keras.optimizers.Adam(), 
              loss=tf.keras.losses.MeanAbsoluteError(),
              metrics=[
                tf.keras.metrics.MeanAbsoluteError(),
                tf.keras.metrics.MeanAbsolutePercentageError(),
              ])

stopping = tf.keras.callbacks.EarlyStopping(
    monitor='val_loss', 
    patience=1000, 
    restore_best_weights=True,
)

model.fit(I, Z, 
    epochs=3000, 
    validation_split=0.3,
    callbacks=[stopping])

model.fit(I, Z, 
    epochs=1, 
    validation_split=0.3,
    callbacks=[stopping])


Z2=model.predict(I)
np.amax(np.abs(Z-Z2))
px.scatter(x=Z.flatten(),y=Z2.flatten()).show()
px.scatter_3d(
    x=np.vstack([V,V]).flatten(),
    y=np.vstack([Tinv,Tinv]).flatten(),
    z=np.vstack([Z,Z2]).flatten(), 
    color=np.vstack([Z*0,Z*0+1]).flatten()).show()

v = V.flatten()
tinv = Tinv.flatten()
z = Z.flatten()
zc = Zc.flatten()


px.scatter_3d(
    x=v,
    y=tinv,
    z=z).show()

mask = np.logical_and(
    np.abs(v-2.0)<0.1, 
    np.abs(zc-0.278)<0.1) #0.244, 0.278, 0.316

px.scatter_3d(
    x=v[mask],
    y=tinv[mask],
    z=z[mask]).show()

px.scatter_3d(
    x=(V).flatten(),
    y=(Tinv).flatten(),
    z=(Z-Z2).flatten()).show()
px.scatter_3d(
    x=(V).flatten(),
    y=(Tinv).flatten(),
    z=(np.abs(Z-Z2)/Z).flatten()).show()
