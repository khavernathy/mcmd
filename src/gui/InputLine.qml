import QtQuick 2.7

Item {
    property string name: "default name"
    property string defaultIn: "default value"
    height: 30
    Text {
       id: thetext
       text: name
       height: 25
       width: 100
    }
    BorderImage {
        id: borderim
        source: "images/lineedit.png"
        border.left: 5; border.top: 5
        border.right: 5; border.bottom: 5
        anchors.left: thetext.right
        width:275
        height:25

    TextInput
    {
        id: editor

        cursorVisible: true;
        font.bold: true
        color: "#151515";
        selectionColor: "Green"
        focus: true
        text: defaultIn

    }
    }
    /*
    TextInput {
       id: theinput
       color: "red"

       anchors.left: thetext.right
       width: 275
       height: 25
       BorderImage {
           id: borderthing
           width: parent.width
           height: parent.height
           source: "images/lineedit.png"
           border.left: 5; border.top: 5
           border.right: 5; border.bottom: 5
       }
       text: defaultIn

    }
    */

}
