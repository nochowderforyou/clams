#ifndef CLAMOURPETITIONMODEL_H
#define CLAMOURPETITIONMODEL_H

#include <QWidget>
#include <QAbstractTableModel>
#include <QStringList>

class CClamour;

class ClamourPetitionModel : public QAbstractTableModel
{
    Q_OBJECT
public:
    explicit ClamourPetitionModel(QWidget *parent = 0);

    enum RowIndex {
        Height = 0,
        TxID = 1,
        PetitionHash = 2,
        URL = 3,
    };

    int rowCount(const QModelIndex &parent) const;
    int columnCount(const QModelIndex &parent) const;
    QVariant data(const QModelIndex &index, int role) const;
    QVariant headerData(int section, Qt::Orientation orientation, int role) const;

    Qt::ItemFlags flags(const QModelIndex &index) const;

    void clear();
    void setPetition(CClamour *newPetition);

signals:

public slots:

private:
    QStringList rowHeaders;
    CClamour *petition;
};

#endif // CLAMOURPETITIONMODEL_H
