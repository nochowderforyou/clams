#ifndef CLAMOURPAGE_H
#define CLAMOURPAGE_H

#include <map>

#include <QWidget>

namespace Ui {
class ClamourPage;
}
class WalletModel;
class CClamour;
class ClamourSupportModel;

QT_BEGIN_NAMESPACE
class QMenu;
QT_END_NAMESPACE

class ClamourPage : public QWidget
{
    Q_OBJECT

public:
    explicit ClamourPage(QWidget *parent = 0);
    ~ClamourPage();

    void setModel(WalletModel *model);

public slots:
    void showClamourTxResult(std::string txID, std::string txError);
    void setClamourSearchResults(CClamour *pResult);
    void showPetitionSupport(std::map<std::string, int> mapSupport);

private slots:
    void on_createPetitionEdit_textChanged();

    void on_createPetitionButton_clicked();

    void on_setVotesButton_clicked();

    void on_searchClamourButton_clicked();

    void on_getPetitionSupportButton_clicked();

    void on_petitionSupportView_customContextMenuRequested(const QPoint &pos);

    void searchHighlightedPetition();

private:
    Ui::ClamourPage *ui;
    WalletModel *model;
    ClamourSupportModel *supportModel;
    QMenu *petitionViewContextMenu;

    void loadVotes();
    void saveVotes();
    void clearSearchTable();
};

#endif // CLAMOURPAGE_H
