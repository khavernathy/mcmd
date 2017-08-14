import QtQuick 2.7

Item {
    property string name: "default name"
    property string defaultIn: "default value"
    Text {
       id: thetext
       text: name
       height: 25
       width: 100
    }
    TextInput {
       id: theinput
       text: defaultIn
       anchors.left: thetext.right
       width: 275
       height: 25
       BorderImage {
           id: borderthing
           source: "images/lineedit.png"
           border.left: 5; border.top: 5
           border.right: 5; border.bottom: 5
       }
    }
    height: 30
}
