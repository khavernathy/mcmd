import QtQuick 2.4

Item {
    width: 1000
    height: 800
    Rectangle {
        id: leftcol
        color: "red"
    }
    Rectangle {
        id: midcol
        anchors.left: leftcol.right
        anchors.right: rightcol.left
        border.color: "black"
    }
    Rectangle {
        id: rightcol
        color: "orange"
    }
}
