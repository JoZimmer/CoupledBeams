from RFEM.enums import  *
from RFEM.initModel import  *
from RFEM.BasicObjects.node import Node
from RFEM.BasicObjects.line import Line
from RFEM.BasicObjects.thickness import Thickness
from RFEM.BasicObjects.surface import Surface
from RFEM.BasicObjects.material import Material
from math import  *


# Schritt 2: Modell initialisieren

Model(new_model=True,model_name="Zahnverbindung")

#SetAddonStatus(Model.clientModel,AddOn.multilayer_surfaces_design_active,True)

Model.clientModel.service.begin_modification("new")

#---------------------------------------------------------------------------------------------------

#Die nachfolgenden Parameter dienen zur Steuerung des Modells"

Zahnlänge =0.5
Zahnhöhe =0.15
Anschnittwinkel =70
Elementdicke =0.16
Modellhöhe =0.8
Holzgüte ="C24"
Längslage = [0,1,0.04]
Querlage = [0,1,0.02]
Lagenaufbau = [Längslage,Querlage,Längslage,Querlage,Längslage]
Lagenwinkel =0

# Knoten aud der Vorderseite

def knoten():
    Node(1,0,0, +Zahnhöhe/2)
    Node(2,0.5*Zahnlänge-(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0,Zahnhöhe/2)
    Node(3,0.5*Zahnlänge+(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0, -Zahnhöhe/2)
    Node(4,1.5*Zahnlänge-(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0, -Zahnhöhe/2)
    Node(5,1.5*Zahnlänge+(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0,Zahnhöhe/2)
    Node(6,2.5*Zahnlänge-(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0,Zahnhöhe/2)
    Node(7,2.5*Zahnlänge+(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0, -Zahnhöhe/2)
    Node(8,3.5*Zahnlänge-(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0, -Zahnhöhe/2)
    Node(9,3.5*Zahnlänge+(Zahnhöhe/2)/tan (radians (Anschnittwinkel)),0, +Zahnhöhe/2)
    Node(10,4*Zahnlänge,0,Zahnhöhe/2)
    Node(11,4*Zahnlänge,0,Modellhöhe/2)
    Node(12,0,0,Modellhöhe/2)
    Node(13,0,0, -Modellhöhe/2)
    Node(14,4*Zahnlänge,0, -Modellhöhe/2)

def linien ():
    Line(1,"1 2")
    Line(2,"2 3")
    Line(3,"3 4")
    Line(4,"4 5")
    Line(5,"5 6")
    Line(6,"6 7")
    Line(7,"7 8")
    Line(8,"8 9")
    Line(9,"9 10")
    Line(19,"10 11")
    Line(20,"11 12")
    Line(21,"12 1")
    Line(22,"1 13")
    Line(23,"13 14")
    Line(24,"10 14")

knoten()
linien()

Material(1, 'C24', params={'material_model': "MODEL_ORTHOTROPIC_2D"})
Thickness.Layers(1, 'BSP', Lagenaufbau)#, params=parameter)

Surface (1,"1 2 3 4 5 6 7 8 9 22 23 24",1)

Model.clientModel.service.finish_modification()
